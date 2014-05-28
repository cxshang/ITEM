
#author: Changxing Shang
#email: cxshang@gmail.com

#The MIT License (MIT)

#Copyright (c) <2013> <Changxing Shang>

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.


from __future__ import print_function
from __future__ import division
import sys
import igraph
cimport cython 
import numpy as np
cimport numpy as np
from time import time

from libc.stdint cimport *
from cityhash cimport *
from fmath cimport *
from popcount cimport *
from libcpp.utility cimport pair
from cython.parallel import prange,threadid
cimport openmp
from arrcontainer cimport *
from stdalgorithm cimport *
from cython.operator cimport dereference as deref, preincrement as inc


DEF FBLEN = 64 #FBLEN is the length of finger bits
DEF USEING_SPARSE_MAPANDSET = False  #indicate whether using google spare map and set or cpp map and set

IF USEING_SPARSE_MAPANDSET:
	from sset cimport sparse_hash_set as sset
	from smap cimport sparse_hash_map as smap
ELSE:
	from libcpp.set cimport set as sset
	from libcpp.map cimport map as smap

ctypedef int64_t MAPKEY_TYPE
ctypedef int32_t MAPVALUE_TYPE
ctypedef int64_t SETKEY_TYPE
ctypedef pair[MAPKEY_TYPE,MAPVALUE_TYPE] MAPAIR
ctypedef smap[MAPKEY_TYPE,MAPVALUE_TYPE] SMAP
ctypedef sset[SETKEY_TYPE] SSET


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(False)
@cython.nonecheck(False)
cdef public void item(char* edgefile, int startindex, float th):
	
	assert(startindex==0 or startindex==1)
	assert(th>=0.0 and th<=1.0)
	
	cdef float esplion=0.001
	cdef int64_t K
	cdef int64_t threadno
	cdef MAPKEY_TYPE key
	cdef MAPAIR kvpair

	ig=igraph.Graph

	cdef int64_t i,j,k,klass
	cdef int64_t maxthreads = openmp.omp_get_max_threads()

	#================================load graph and get csr format of graph=========================================
	print( 'loading %s ...' % edgefile, file=sys.stderr)
	time0=time()
	
	g=ig.Read_Edgelist(edgefile ,directed=False)	

	cdef int64_t nodescnt=g.vcount()
	cdef int64_t edgescnt=g.ecount()
	
	xtemp=[edge.tuple[0] for edge in g.es]
	cdef int64_t[::1] pts1=np.array(xtemp)
	xtemp=[edge.tuple[1] for edge in g.es]
	cdef int64_t[::1] pts2=np.array(xtemp)

	cdef int64_t eidx
	cdef int64_t tempswitch
	for eidx in prange(edgescnt ,nogil=True, num_threads=maxthreads):
		if pts1[eidx]>pts2[eidx]:
			tempswitch=pts1[eidx]
			pts1[eidx]=pts2[eidx]
			pts2[eidx]=tempswitch
	
	cdef int64_t[::1]  degrees1= np.take(g.degree(), pts1)
	cdef int64_t[::1]  degrees2= np.take(g.degree(), pts2)
	cdef int64_t[::1]  degrees= np.take(g.degree(), range(nodescnt))
	cdef float avgdegree=sum(degrees)/nodescnt

	cdef int64_t  cnts=sum(np.minimum(degrees1,degrees2)+2)
	cdef int64_t[::1]  edgesidx=np.zeros(edgescnt+1, np.int64, order='C') 
	cdef int64_t[::1]  nodesidx=np.zeros(cnts, np.int64, order='C')
	cdef float[::1] values=np.zeros(cnts, np.float32, order='C')

	cdef float[::1] idfs=np.zeros(nodescnt, np.float32, order='C') #term's invert doc frequency
	
	for i in prange(nodescnt ,nogil=True, num_threads=maxthreads):
		if degrees[i]!=0:
			idfs[i]=log(edgescnt)-log(degrees[i])
			
	cdef int64_t[::1] dfs=np.zeros(nodescnt, np.int64, order='C') #term doc frequency
	
	cdef int64_t curpos
	cdef int64_t ndidx,	nbndidx
	cdef int64_t idx1
	cdef int64_t idx2
	cdef int64_t startidx
	cdef int64_t endidx
	cdef int8_t pts1inserted
	cdef int8_t pts2inserted
	
	graphadjlist=g.get_adjlist()
	cdef int64_t[::1]  adjlist=np.zeros(np.sum(degrees), np.int64, order='C')
	cdef int64_t[::1]  adjlistdivs=np.zeros(nodescnt+1, np.int64, order='C')
	j=0
	for i in range(nodescnt):
		for ndidx in graphadjlist[i]:
			adjlist[j]=ndidx
			j+=1
		adjlistdivs[i+1]=j
	
	curpos=0
	for eidx in range(edgescnt):
		
		idx1=adjlistdivs[pts1[eidx]]
		idx2=adjlistdivs[pts2[eidx]]
		pts1inserted=0
		pts2inserted=0
		while ((idx1<adjlistdivs[pts1[eidx]+1]) and (idx2<adjlistdivs[pts2[eidx]+1])):
			if adjlist[idx1]==adjlist[idx2]:
				nodesidx[curpos]=adjlist[idx1]
				values[curpos]=idfs[adjlist[idx1]]
				dfs[adjlist[idx1]]+=1
				curpos+=1
				idx1+=1
				idx2+=1
			elif adjlist[idx2]==pts1[eidx] and adjlist[idx2]<adjlist[idx1] :
				nodesidx[curpos]=pts1[eidx]
				values[curpos]=idfs[pts1[eidx]]
				dfs[pts1[eidx]]+=1
				curpos+=1
				idx2+=1
				pts1inserted=1
			elif adjlist[idx1]==pts2[eidx] and adjlist[idx1]<adjlist[idx2] :
				nodesidx[curpos]=pts2[eidx]
				values[curpos]=idfs[pts2[eidx]]
				dfs[pts2[eidx]]+=1
				curpos+=1
				idx1+=1			
				pts2inserted=1
			elif adjlist[idx1]<adjlist[idx2]:
				idx1+=1
			elif adjlist[idx1]>adjlist[idx2]:
				idx2+=1
				
		if pts1inserted==0:
			nodesidx[curpos]=pts1[eidx]
			values[curpos]=idfs[pts1[eidx]]
			dfs[pts1[eidx]]+=1
			curpos+=1
		if pts2inserted==0:
			nodesidx[curpos]=pts2[eidx]
			values[curpos]=idfs[pts2[eidx]]
			dfs[pts2[eidx]]+=1
			curpos+=1				
		edgesidx[eidx+1]=curpos
		
	assert(curpos<cnts)
	
	print( "edges count is %d"%edgescnt, file=sys.stderr)
	print( "nodes count is %d"% nodescnt, file=sys.stderr)
	print( "average degree per node is %.2f"% (sum(degrees)/nodescnt), file=sys.stderr)
	print( "average shared nodes (plus two endpoints) per edge is %.2f"% (curpos/edgescnt) , file=sys.stderr)
	print( "sparse ratio of jaccard matrix is %.5f"%(curpos/(edgescnt*nodescnt)), file=sys.stderr)
	print( "max threads count is %d"%maxthreads, file=sys.stderr)
	
	#jaccardmtx=sparse.csr_matrix( (np.ones(curpos, np.uint32,order='C'),nodesidx[:curpos],edgesidx), shape=(edgescnt,nodescnt) )
		
	#=====================================cacluate fingerprint of each edge====================================
	time1=time()

	IF FBLEN==32:
		cdef uint32_t[::1] edgefingers=np.zeros(edgescnt, np.uint32, order='C')
		cdef uint32_t hashval
	IF FBLEN==64:
		cdef uint64_t[::1] edgefingers=np.zeros(edgescnt, np.uint64, order='C')
		cdef uint64_t hashval
		
	cdef int64_t nnzpt

	cdef float[::1] simhash=np.zeros(FBLEN, np.float32, order='C') 
	cdef float negsum
	for eidx in range(edgescnt):
		negsum=0.0
		startidx=edgesidx[eidx]
		endidx=edgesidx[eidx+1]
		
		for j in range(startidx, endidx):
			negsum=negsum-values[j]
		for j in range(FBLEN):
			simhash[j]=negsum

		for i in range(startidx, endidx): 
			nnzpt=nodesidx[i]
			IF FBLEN==32:
				hashval=cityhash32(nnzpt) 
			IF FBLEN==64:
				hashval=cityhash64(nnzpt) 	

			for k in range(FBLEN): 
				if (hashval>>k)&1:
					simhash[k]=simhash[k]+(values[i]+values[i])
					
		for k in range(FBLEN): 
			if (simhash[k]>0):
				edgefingers[eidx]|=(1<<k)

	print( "Cost %.2f ms for finger." % ((time()-time1)*1000), file=sys.stderr)

	#====================computing strength of each edge==============================
	time2=time()
	cdef int64_t maxdegree
	cdef float[::1] strength=np.empty(edgescnt, np.float32, order='C') 
	for eidx in prange(edgescnt ,nogil=True, num_threads=maxthreads):
		maxdegree=max(degrees1[eidx], degrees2[eidx])
		strength[eidx]=(edgesidx[eidx+1]-edgesidx[eidx]-1)/maxdegree
	print( "Cost %.2f us for strength." % ((time()-time2)*1000000), file=sys.stderr)
	
	#====================decide whether a node is zero order==============================
	time3=time()
	graphinclist=g.get_inclist()
	cdef int64_t[::1] inclist=np.zeros(np.sum(degrees), np.int64, order='C')
	cdef int64_t[::1] inclistdivs=np.zeros(nodescnt+1, np.int64, order='C')
	j=0
	for ndidx in range(nodescnt):
		for eidx in graphinclist[ndidx]:
			inclist[j]=eidx
			j+=1
		inclistdivs[ndidx+1]=j
	print( "Cost %.2f ms for init inclist." % ((time()-time3)*1000), file=sys.stderr)
	
	cdef int64_t[::1] zerordernodes=np.zeros(nodescnt, np.int64, order='C')
	for ndidx in range(nodescnt):
		for i in range(inclistdivs[ndidx], inclistdivs[ndidx+1]):
			eidx=inclist[i]
			zerordernodes[ndidx]+=edgesidx[eidx+1]-edgesidx[eidx]-2
	
	#====================computing nodespecificity of each nodes==============================
	time4=time()
	cdef float[::1] nodespecificity=np.zeros(nodescnt, np.float32, order='C') 
	cdef float sumofstrength
	for ndidx in prange(nodescnt ,nogil=True, num_threads=maxthreads):
		sumofstrength=0.0
		for i in range(inclistdivs[ndidx], inclistdivs[ndidx+1]) :
			eidx=inclist[i]
			sumofstrength=sumofstrength+strength[eidx]
		nodespecificity[ndidx]=sumofstrength/(inclistdivs[ndidx+1]-inclistdivs[ndidx])
	print( "Cost %.2f us for nodespecificity." % ((time()-time4)*1000000), file=sys.stderr)
	
	#====================computing edgespecificity of each edges==============================
	time5=time()
	cdef float[::1] edgespecificity=np.empty(edgescnt, np.float32, order='C') 
	for eidx in prange(edgescnt ,nogil=True, num_threads=maxthreads):
		edgespecificity[eidx]=min(nodespecificity[pts1[eidx]], nodespecificity[pts2[eidx]])
	print( "Cost %.2f us for edgespecificity." % ((time()-time5)*1000000), file=sys.stderr)
	
	#====================computing reputation of each edge, reputation is sum of similarity==============================
	time6=time()
	cdef int64_t nbedgeidx
	cdef float[::1] reputation=np.zeros(edgescnt, np.float32, order='C')
	for eidx in prange(edgescnt ,nogil=True, num_threads=maxthreads):
		reputation[eidx]=-2*FBLEN
		for i in range(inclistdivs[pts1[eidx]],inclistdivs[pts1[eidx]+1]):
			nbedgeidx=inclist[i]
			reputation[eidx]+=(FBLEN-hammdist(edgefingers[eidx], edgefingers[nbedgeidx]))
		for i in range(inclistdivs[pts2[eidx]],inclistdivs[pts2[eidx]+1]):
			nbedgeidx=inclist[i]
			reputation[eidx]+=(FBLEN-hammdist(edgefingers[eidx], edgefingers[nbedgeidx]))
	print( "Cost %.2f ms for reputation." % ((time()-time6)*1000), file=sys.stderr)
	
	
	#================cacluate isclimaxs of each edge===============================
	time7=time()
	cdef float[::1] altitude=np.zeros(edgescnt, np.float32, order='C') 
	for eidx in range(edgescnt):
		altitude[eidx]=reputation[eidx]*strength[eidx]*edgespecificity[eidx]
	cdef int8_t[::1] isclimaxs=np.zeros(edgescnt, np.int8, order='C') 
	for eidx in prange(edgescnt ,nogil=True, num_threads=maxthreads):
		if edgesidx[eidx+1]-edgesidx[eidx]>=4: 
			isclimaxs[eidx]=1
	for eidx in range(edgescnt):
		if isclimaxs[eidx]==0:
			continue
		for i in range(inclistdivs[pts1[eidx]],inclistdivs[pts1[eidx]+1]):
			nbedgeidx=inclist[i]
			if altitude[eidx]<altitude[nbedgeidx]:
				isclimaxs[eidx]=0
				break
			elif altitude[eidx]>altitude[nbedgeidx]:
				isclimaxs[nbedgeidx]=0
			elif eidx<nbedgeidx:
				isclimaxs[nbedgeidx]=0
			elif eidx>nbedgeidx:
				isclimaxs[eidx]=0
				
		if isclimaxs[eidx]==0:
			continue
		for i in range(inclistdivs[pts2[eidx]],inclistdivs[pts2[eidx]+1]):
			nbedgeidx=inclist[i]
			if altitude[eidx]<altitude[nbedgeidx]:
				isclimaxs[eidx]=0
				break
			elif altitude[eidx]>altitude[nbedgeidx]:
				isclimaxs[nbedgeidx]=0
			elif eidx<nbedgeidx:
				isclimaxs[nbedgeidx]=0
			elif eidx>nbedgeidx:
				isclimaxs[eidx]=0
				
	print( "Cost %.2f ms for climaxs." % ((time()-time7)*1000), file=sys.stderr)
									
	
	cdef int64_t[::1] candidatedges=np.where(isclimaxs)[0]
	cdef int64_t candidategdescnt=len(candidatedges)
	print( "at first, total class(climaxs)=%d" %candidategdescnt, file=sys.stderr)
	
	#cdef np.ndarray[int8_t,ndim=1, mode="c"] ishubnode=np.zeros(nodescnt, np.int8, order='C') 
	#argsortofspecificity=np.argsort(nodespecificity)
	#for i in range(int(nodescnt*0.05)):
	#	ishubnode[argsortofspecificity[i]]|=1
		
	#=====================================selecting  seeds ====================================

	cdef int8_t[::1] selected=np.zeros(candidategdescnt, np.int8, order='C') 	
	cdef int8_t[::1] selectedandmarked=np.zeros(candidategdescnt, np.int8, order='C') 	
	cdef int64_t numofalonecommittees=0
	
	cdef int64_t[::1] committeenodes=np.zeros(edgescnt*2, np.int64, order='C')
	cdef int64_t[::1] committeedivs=np.zeros(candidategdescnt+1, np.int64, order='C')	
	cdef int64_t[::1] committeeweits=np.zeros(candidategdescnt, np.int64, order='C')
	cdef float[::1] entropyofcommittees=np.zeros(candidategdescnt, np.float32, order='C') 

	curpos=0
	for k in range(candidategdescnt):
		if selectedandmarked[k]:
			committeedivs[k]=curpos
			committeedivs[k+1]=curpos
			continue
		
		eidx=candidatedges[k]
		startidx=edgesidx[eidx]
		endidx=edgesidx[eidx+1]
		committeedivs[k]=curpos
		for i in range(startidx, endidx):
			ndidx=nodesidx[i]
			committeenodes[curpos]=ndidx
			curpos+=1
		committeedivs[k+1]=curpos
			
		committeeweits[k]=endidx-startidx
		entropyofcommittees[k]=log(committeeweits[k])
	
	cdef arrcontainer[SSET] nodes2communitieset=deref(new arrcontainer[SSET](nodescnt))
	cdef SMAP communitiesoverlapnodestj
	for k in range(candidategdescnt):
		if selectedandmarked[k]:
			continue	
		for i in range(committeedivs[k],committeedivs[k+1]):
			ndidx=committeenodes[i]
			nodes2communitieset[ndidx].insert(k)
	
	time11=time()
	cdef int64_t[::1] predistr=np.zeros(nodescnt, np.int64, order='C')
	cdef float[::1] entropyofpredistr=np.zeros(nodescnt, np.float32, order='C')
	cdef float sumoflogpredistr=0
	cdef int64_t presumofweits=0
	cdef int64_t pickededgeidx
	
	#++++++++get first seed++++++++++++
	pickededgeidx=np.argmax(entropyofcommittees) 
	selected[pickededgeidx]=1
	selectedandmarked[pickededgeidx]=1
	startidx=committeedivs[pickededgeidx]
	endidx=committeedivs[pickededgeidx+1]
	for i in range(startidx, endidx):
		predistr[committeenodes[i]]+=1
		entropyofpredistr[committeenodes[i]]=-predistr[committeenodes[i]]*log(predistr[committeenodes[i]])
		sumoflogpredistr+=entropyofpredistr[committeenodes[i]]
		presumofweits+=1
	
	#++++++++get last K-1 seeds++++++++++++
	cdef float[::1] gigdifference=np.empty(candidategdescnt, np.float32, order='C')
	cdef int64_t cntselectedandmarked=1+numofalonecommittees
	cdef int64_t cntglobaloverlapnodes
	while cntselectedandmarked<candidategdescnt:

		communitiesoverlapnodestj.swap(deref(new SMAP()))
		for i in range(committeedivs[pickededgeidx], committeedivs[pickededgeidx+1]):
			ndidx=committeenodes[i]
			for k in nodes2communitieset[ndidx]:
				communitiesoverlapnodestj[k]+=1
		for kvpair in communitiesoverlapnodestj:
			k=kvpair.first
			if selectedandmarked[k]:
				continue
			cntoverlapnodes=kvpair.second
			if cntoverlapnodes/(committeedivs[k+1]-committeedivs[k])>th:
				selectedandmarked[k]=1
				cntselectedandmarked+=1
		
		for k in prange(candidategdescnt ,nogil=True, num_threads=maxthreads): 
			if selectedandmarked[k]:
				continue
			
			startidx=committeedivs[k]
			endidx=committeedivs[k+1]
			
			gigdifference[k]=sumoflogpredistr-entropyofcommittees[k]
			for j in range(startidx, endidx):
				ndidx=committeenodes[j]
				gigdifference[k]-=entropyofpredistr[ndidx]
				gigdifference[k]-=(predistr[ndidx]+1)*log(predistr[ndidx]+1)
			gigdifference[k]+=(committeeweits[k]+presumofweits)*log(committeeweits[k]+presumofweits)
		
		maxvalues=-2**10
		pickededgeidx=-1
		for k in range(candidategdescnt):
			if selectedandmarked[k]:
				continue
			elif gigdifference[k]>maxvalues:
				maxvalues=gigdifference[k]
				pickededgeidx=k
		
		if pickededgeidx>=0:
			selected[pickededgeidx]=1
			selectedandmarked[pickededgeidx]=1
			cntselectedandmarked+=1
			startidx=committeedivs[pickededgeidx]
			endidx=committeedivs[pickededgeidx+1]
			for j in range(startidx, endidx):
				predistr[committeenodes[j]]+=1
				presumofweits+=1
				
				sumoflogpredistr-=entropyofpredistr[committeenodes[j]]
				entropyofpredistr[committeenodes[j]]=-predistr[committeenodes[j]]*log(predistr[committeenodes[j]])
				sumoflogpredistr+=entropyofpredistr[committeenodes[j]]

	K=0
	for k in range(candidategdescnt):
		K+=selected[k]
	print( "Cost %.2f second for select seed MGIG." % ((time()-time11)), file=sys.stderr)
	print("after filtering, there are %d edge groups left as seeds."%K, file=sys.stderr)		
	
	#=====================construct committees==========================
	time9=time()
	cdef arrcontainer[SSET] committees=deref(new arrcontainer[SSET](candidategdescnt))
	
	for k in range(candidategdescnt):
		if not selected[k]:
			continue
		
		eidx=candidatedges[k]
		startidx=edgesidx[eidx]
		endidx=edgesidx[eidx+1]
		
		for i in range(startidx, endidx):
			for j in range(i+1, endidx):
				idx1=g.get_eid(nodesidx[i],nodesidx[j], directed=False, error=False)
				if idx1 != -1 :
					committees[k].insert(idx1) 
					
	print( "Cost %.2f ms for making committees." % ((time()-time9)*1000), file=sys.stderr)
	
	#==========================generate Prob(node|edge) in values=======================
	cdef int64_t maxdfs=np.max(dfs)
	for j in prange(nodescnt ,nogil=True, num_threads=maxthreads):
		idfs[j]=log(maxdfs+1) - log(dfs[j])
	
	for i in prange(edgescnt ,nogil=True, num_threads=maxthreads):
		startidx=edgesidx[i]
		endidx=edgesidx[i+1]
		for j in range(startidx, endidx):
			values[j]=dfs[nodesidx[j]]*idfs[nodesidx[j]]

	#======================assign labels to init seed edges=============================
	cdef int64_t[::1] labels=np.empty(edgescnt, np.int64, order='C')
	for i in prange(edgescnt ,nogil=True, num_threads=maxthreads):
		labels[i]=-1
	cdef int64_t classlabel=0
	for j in range(candidategdescnt):
		if selected[j]==0:
			continue
		for eidx in committees[j]:
			labels[eidx]=classlabel
		classlabel+=1
	
	#====================================================================================
	cdef int64_t[::1] centroidrows=np.zeros(K+1, np.int64, order='C')
	cdef np.ndarray[np.int64_t,ndim=1, mode="c"] npcentroidcols=np.zeros(edgesidx[edgescnt], np.int64, order='C')
	cdef int64_t[::1] centroidcols=npcentroidcols
	cdef float[::1] logcentroidvalues=np.zeros(edgesidx[edgescnt], np.float32, order='C')
	cdef int64_t[::1] communitysizes=np.zeros(K, np.int64, order='C')
	cdef float[::1] kldivergences=np.empty(edgescnt, np.float32, order='C')
	
	cdef int64_t iteridx
	cdef int64_t startidxklass
	cdef int64_t endidxklass
	cdef int8_t[::1] inclusivenodes=np.zeros(nodescnt, np.int8, order='C')

	cdef float[::1] sumprob=np.zeros(K, np.float32, order='C') 
	cdef SMAP bigmap
	cdef arrcontainer[SMAP] arrsmap=deref(new arrcontainer[SMAP](maxthreads))
	cdef SMAP nodesoccurs
		
	cdef int64_t cntcurcommunity
	cdef np.float32_t weights
	cdef np.float32_t logsmoothing
	cdef int64_t cntchangededges, oldcntchangededges=0
	cdef np.float32_t kldiv
	cdef np.float32_t unmatchedweit
	cdef int64_t totalcommunitysizes
	cdef int64_t nnzentries
	cdef int64_t multiplier
	cdef int64_t cntsharenodes
	cdef int64_t innerpoint
	cdef int64_t outerpoint
	cdef int64_t nodemask
	cdef int64_t pre_k
	cdef np.float32_t priorprobclass
	cdef int64_t cntintersectionnodes
	
	cdef int64_t adjedgeidx
	cdef int64_t adjedgecnt
	cdef int64_t[::1] adjedgesofclass=np.zeros(edgescnt, np.int64, order='C')
	cdef SSET setofadjedgesofclass
	
	cdef int8_t eleftrorbits=int(log(nodescnt)/log(2))
	if (1<<eleftrorbits)<nodescnt:
		eleftrorbits+=1
	if (1<<eleftrorbits)<nodescnt:
		eleftrorbits+=1
	if (1<<eleftrorbits)<nodescnt:
		eleftrorbits+=1
	assert((1<<eleftrorbits)>=nodescnt)
	assert((((K-1)<<eleftrorbits)|(nodescnt-1)) < ((1<<63)))
	nodemask=(1<<eleftrorbits)-1

	cdef int64_t blocksize
	if edgescnt%maxthreads==0:
		blocksize=edgescnt//maxthreads
	else:
		blocksize=edgescnt//maxthreads+1	

	cdef int64_t mididx, leftidx, rightidx, newstartidxklass
	#=====================================EM algorithm ===========================================
	print( 'start EM iters...', file=sys.stderr)

	iteridx=0
	for j in prange(edgescnt, nogil=True, num_threads=maxthreads):
		kldivergences[j]=(2**21)*1.0
	
	while True:
		iteridx+=1
		
		#++++++++++++++++++++++++++++++++Maximization step+++++++++++++++++++++++++
		bigmap.swap(deref(new SMAP()))
		for threadno in range(maxthreads):
			arrsmap[threadno].swap(deref(new SMAP()))
		for k in range(K):
			sumprob[k]=0.0
		for klass in range(K):
			communitysizes[klass]=0
		totalcommunitysizes=0
		for k in range(K): 
			centroidrows[k]=-1

		for eidx in range(edgescnt):
			if labels[eidx] >=0: 
				communitysizes[labels[eidx]]+=1
				totalcommunitysizes+=1
		
		for eidx in prange(edgescnt, nogil=True, num_threads=maxthreads):
			if labels[eidx] >=0: 
				threadno=eidx//blocksize
				startidx=edgesidx[eidx]
				endidx=edgesidx[eidx+1]
				for j in range(startidx, endidx):
					key=(labels[eidx]<<eleftrorbits)|nodesidx[j]
					arrsmap[threadno][key]+=1
		
		nodesoccurs.swap(deref(new SMAP()))
		for threadno in range(maxthreads):
			if arrsmap[threadno].empty():
				continue
			for kvpair in arrsmap[threadno]:
				key=kvpair.first
				bigmap[key]+=arrsmap[threadno][key]
				
				ndidx=(key & nodemask)
				nodesoccurs[ndidx]+=arrsmap[threadno][key]
				
		#assert(bigmap.size()<edgesidx[edgescnt])
		
		
		nnzentries=0
		for kvpair in bigmap:
			centroidcols[nnzentries]=kvpair.first
			nnzentries+=1
		npcentroidcols[0:nnzentries].sort()
		#assert(nnzentries<=bigmap.max_size())	
		
		pre_k=-1
		cntcurcommunity=0
		for j in range(nnzentries):
			key=centroidcols[j]
			k=(key>>eleftrorbits)
			if k != pre_k:
				centroidrows[k]=j
				cntcurcommunity+=1
			pre_k=k
		centroidrows[K]=nnzentries
		for k in range(K-1,-1,-1): 
			if centroidrows[k]==-1:
				centroidrows[k]=centroidrows[k+1]	
		
		for j in range(nnzentries):
			key=centroidcols[j]
			ndidx=(key & nodemask)
			k=(key>>eleftrorbits)
			centroidcols[j]=ndidx
			logcentroidvalues[j]=(bigmap[key]/communitysizes[k])/((nodesoccurs[ndidx]-bigmap[key]+1)/(totalcommunitysizes-communitysizes[k]))
			sumprob[k]+=logcentroidvalues[j]

		logsmoothing=0.0
		for k in range(K):
			if communitysizes[k]==0:
				continue
			startidx=centroidrows[k]
			endidx=centroidrows[k+1]
			sumprob[k]=log(sumprob[k])
			for j in range(startidx, endidx):
				logcentroidvalues[j]=-( log(logcentroidvalues[j])-sumprob[k]) # -log( Porb(node|cluster) )
				if logsmoothing<logcentroidvalues[j]:
					logsmoothing=logcentroidvalues[j]
		logsmoothing=logsmoothing*1.001
		#assert(logsmoothing<2**20)
		
		for j in prange(edgescnt, nogil=True, num_threads=maxthreads):
			if labels[j]<0:
				continue
				
			startidxklass=centroidrows[labels[j]]
			endidxklass=centroidrows[labels[j]+1]
			startidx=edgesidx[j]
			endidx=edgesidx[j+1]
			
			newstartidxklass=startidxklass
			unmatchedweit=0
			kldiv=-( log(communitysizes[labels[j]]) - log(totalcommunitysizes) )
			for i in range(startidx, endidx):
				leftidx=newstartidxklass
				rightidx=endidxklass-1
				
				while (leftidx<=rightidx):
					mididx=(leftidx+rightidx)//2
					if nodesidx[i]==centroidcols[mididx]:
						kldiv=kldiv+values[i]*logcentroidvalues[mididx]
						newstartidxklass=mididx+1
						break
					elif nodesidx[i]<centroidcols[mididx]:
						rightidx=mididx-1
					else:
						leftidx=mididx+1
				else:
					unmatchedweit=unmatchedweit+values[i]
			
			kldivergences[j]=kldiv+unmatchedweit*logsmoothing
			
		#++++++++++++++++++++++++++++++++Expectation step +++++++++++++++++++++++++
		cntchangededges=0
		for klass in range(K):
			if communitysizes[k]==0:
				continue

			for j in prange(nodescnt, nogil=True, num_threads=maxthreads):
				inclusivenodes[j]=0
			for j in range(edgescnt):
				if labels[j]==klass:
					inclusivenodes[pts1[j]]|=1
					inclusivenodes[pts2[j]]|=1
			
			setofadjedgesofclass.swap(deref(new SSET()))
			for ndidx in range(nodescnt):
				if inclusivenodes[ndidx]:
					for i in range(inclistdivs[ndidx],inclistdivs[ndidx+1]):
						eidx=inclist[i]
						setofadjedgesofclass.insert(eidx)
			adjedgecnt=0
			for eidx in setofadjedgesofclass:
				adjedgesofclass[adjedgecnt]=eidx
				adjedgecnt+=1
					
			startidxklass=centroidrows[klass]
			endidxklass=centroidrows[klass+1]
			
			priorprobclass=-( log(communitysizes[klass]) - log(totalcommunitysizes) )
			for adjedgeidx in prange(adjedgecnt, nogil=True, num_threads=maxthreads):
				eidx=adjedgesofclass[adjedgeidx]
				
				if klass==labels[eidx]:
					continue
				
				startidx=edgesidx[eidx]
				endidx=edgesidx[eidx+1]
				
				if not (inclusivenodes[pts1[eidx]] and inclusivenodes[pts2[eidx]]):
					if inclusivenodes[pts1[eidx]] :
						innerpoint=pts1[eidx]
						outerpoint=pts2[eidx]
					else:
						innerpoint=pts2[eidx]
						outerpoint=pts1[eidx]
						
					cntintersectionnodes=0
					for i in range(startidx, endidx):
						cntintersectionnodes=cntintersectionnodes+inclusivenodes[nodesidx[i]]
					if endidx-startidx>2 and cntintersectionnodes<=1:
						continue
					elif endidx-startidx==2:	
						if zerordernodes[outerpoint]>0:
							continue
						elif (cntintersectionnodes<<1)<degrees[outerpoint]:
							continue
					
				unmatchedweit=0
				kldiv=priorprobclass
				newstartidxklass=startidxklass
				for i in range(startidx, endidx):
					leftidx=newstartidxklass
					rightidx=endidxklass-1
					while (leftidx<=rightidx):
						mididx=(leftidx+rightidx)//2
						if nodesidx[i]==centroidcols[mididx]:
							kldiv=kldiv+values[i]*logcentroidvalues[mididx]
							newstartidxklass=mididx+1
							break
						elif nodesidx[i]<centroidcols[mididx]:
							rightidx=mididx-1
						else:
							leftidx=mididx+1
					else:
						unmatchedweit=unmatchedweit+values[i]
						
					if kldiv>kldivergences[eidx] : 
						break 
					
				kldiv=kldiv+unmatchedweit*logsmoothing
				if kldiv<kldivergences[eidx] : 
					cntchangededges+=1
					kldivergences[eidx]=kldiv
					labels[eidx]=klass
		
		print( 'there are %d edges labels changed in iter %d'%(cntchangededges, iteridx), file=sys.stderr)	
		if cntchangededges<max(edgescnt*esplion, 3) or iteridx>=20 or \
		abs(cntchangededges-oldcntchangededges)<0.05*cntchangededges :
			break
		else:
			oldcntchangededges=cntchangededges
		
	#+++++++++++++++++++++++++++processing unlabeled nodes+++++++++++++++++++++++++++++++++++	
	print( 'EM iters finished', file=sys.stderr)

	bigmap.swap(deref(new SMAP()))
	for threadno in range(maxthreads):
		arrsmap[threadno].swap(deref(new SMAP()))

	cdef int64_t unlabelededges=0
	cdef arrcontainer[SSET] nodescolors=deref(new arrcontainer[SSET](nodescnt))
	for i in range(edgescnt):
		if labels[i]>=0:
			nodescolors[pts1[i]].insert(labels[i])
			nodescolors[pts2[i]].insert(labels[i])
		else:
			unlabelededges+=1
	print( "there are %d unlabeled edges."%unlabelededges , file=sys.stderr)		
	
	cdef int64_t cntunlabelednodes=0
	for ndidx in range(nodescnt):
		if nodescolors[ndidx].size()==0:
			cntunlabelednodes+=1
	print( "there are %d unlabeled nodes."%cntunlabelednodes , file=sys.stderr)
	
	
	#==========================processing unlabeled nodes=====================================
	print( "start processing unlabeled nodes...", file=sys.stderr)
	cdef arrcontainer[SMAP] cntcolorsdiffersmaps=deref(new arrcontainer[SMAP](nodescnt)) 
	cdef arrcontainer[SMAP] cntneibercolorsmaps=deref(new arrcontainer[SMAP](nodescnt)) 
	cdef int64_t[::1] maxcntofneibercolors=np.zeros(nodescnt, np.int64, order='C')
	IF USEING_SPARSE_MAPANDSET:
		for ndidx in range(nodescnt):
			cntcolorsdiffersmaps[ndidx].set_deleted_key(-1)
	
	for ndidx in range(nodescnt):
		for i in range(adjlistdivs[ndidx], adjlistdivs[ndidx+1]):
			nbndidx=adjlist[i]
			for k in nodescolors[nbndidx]:
				cntneibercolorsmaps[ndidx][k]+=1
				if maxcntofneibercolors[ndidx]<cntneibercolorsmaps[ndidx][k]:
					maxcntofneibercolors[ndidx]=cntneibercolorsmaps[ndidx][k]
				if nodescolors[ndidx].find(k)==nodescolors[ndidx].end():
					cntcolorsdiffersmaps[ndidx][k]+=1

	iteridx=0	
	cdef np.int64_t changed=1
	cdef int64_t[::1] newaddcolors=np.empty(nodescnt, np.int64, order='C')
	while changed :
		iteridx+=1
		changed=0
			
		for ndidx in range(nodescnt):
			newaddcolors[ndidx]=-1
		#add new color according colors of neighbors
		for ndidx in range(nodescnt):
			for kvpair in cntcolorsdiffersmaps[ndidx]:
				if kvpair.second==maxcntofneibercolors[ndidx]:
					changed+=1
					newaddcolors[ndidx]=kvpair.first
					break
					
		print( "there are %d nodes changed labels in iter %d."%(changed, iteridx) , file=sys.stderr)
		if changed==0:
			break
			
		for ndidx in range(nodescnt):	
			k=newaddcolors[ndidx]
			if k==-1:
				continue
			nodescolors[ndidx].insert(k)
			cntneibercolorsmaps[ndidx][k]=cntcolorsdiffersmaps[ndidx][k]
			if maxcntofneibercolors[ndidx]<cntneibercolorsmaps[ndidx][k]:
				maxcntofneibercolors[ndidx]=cntneibercolorsmaps[ndidx][k]
			cntcolorsdiffersmaps[ndidx].erase(k)	
				
		#change colors impact neighbors 
		for ndidx in range(nodescnt):
			k=newaddcolors[ndidx]
			if k==-1:
				continue
			for i in range(adjlistdivs[ndidx], adjlistdivs[ndidx+1]):
				nbndidx=adjlist[i]
				cntneibercolorsmaps[nbndidx][k]+=1
				if maxcntofneibercolors[nbndidx]<cntneibercolorsmaps[nbndidx][k]:
					maxcntofneibercolors[nbndidx]=cntneibercolorsmaps[nbndidx][k]
				if nodescolors[nbndidx].find(k)==nodescolors[nbndidx].end():
					cntcolorsdiffersmaps[nbndidx][k]+=1
	
	#==========================output results=====================================
	cntunlabelednodes=0
	for ndidx in range(nodescnt):
		if nodescolors[ndidx].size()==0:
			cntunlabelednodes+=1
	print( "finally, there are %d unlabeled nodes."%cntunlabelednodes , file=sys.stderr)
		
	cdef arrcontainer[SSET] nodespercommunity=deref(new arrcontainer[SSET](K))
	for ndidx in range(nodescnt):
		for k in nodescolors[ndidx]:
			nodespercommunity[k].insert(ndidx) 
	klass=K
	for k in range(K):
		if nodespercommunity[k].size()==0:
			klass-=1
			continue
		if startindex==0:
			print(  "%s" % ' '.join([str(ndidx) for ndidx in nodespercommunity[k]])  )
		elif startindex==1:
			print(  "%s" % ' '.join([str(ndidx+1) for ndidx in nodespercommunity[k]])  ) 
			
	print( "there are total %d communities."%klass , file=sys.stderr)
	print( "totally cost %.2f seconds. Finished." % ((time()-time0)), file=sys.stderr)
	
	
