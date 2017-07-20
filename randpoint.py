import numpy as np
from math import *
import matplotlib.pyplot as plt
import time
import sys
from itertools import combinations
import random

XMARGIN=1.5
EPS=0.00001
POINT_SIZE=0.05
N=10
Q=6


#	
# construct a random sample of n axis-parallel squares
#
def randpoints(n,k,seed):  
	if seed != -1:
	    np.random.seed(seed)

	xhi=1/2*(1+k)/(sqrt(1+k**2))
	xlow=-xhi
	yhi=np.sqrt(2)/4*N
	ylow=-yhi

	x=np.random.uniform(xlow,xhi,n)
	y=np.sort(np.random.uniform(ylow,yhi,n))
	
	return np.matrix([[x[i],y[i]] for i in range(n)])


# #
# #  draw axis-parallel squares
# #  with a given centers
# #
def draw_squares(n,k,centers,rotator, P):
		
	xhi=0.5*(1+k)/sqrt(1+k**2)
	xlow=-xhi

	yhi=sqrt(2)/4*N
	ylow=-yhi

	Diag=np.matrix([[0,yhi],[0,ylow]]).transpose()
	Diag_transp=(rotator.dot(Diag)).transpose()
	d_line=plt.Polygon(Diag_transp,closed=False,fill=False, color='g')
	

	plt.axes()
	plt.gca().add_patch(d_line)

	# draw the squares
	for i in range(n):
		LA=centers-[0.5,0.5]
		la=np.ravel(LA[i,:])
		sq=plt.Rectangle((la[0],la[1]),1,1,fill=False, color='b')
		plt.gca().add_patch(sq)

	# draw orthogonal lines	
	Orth=[[[-XMARGIN,ylow+sqrt(2)/2*i],[XMARGIN,ylow+sqrt(2)/2*i]] for i in range(N+1)]
	for i in range(N+1):
		ort=np.matrix(Orth[i]).transpose()
		ort_transp=(rotator.dot(ort)).transpose()
		ort_line=plt.Polygon(ort_transp,closed=False, linestyle='dashed')
		plt.gca().add_patch(ort_line)

	for i in range(len(P)):
		p=plt.Circle((P[i][0],P[i][1]),POINT_SIZE,color='r')
		plt.gca().add_patch(p)

	plt.axis('scaled')
	plt.show()

# def make_ground_set(square_centers):

# def square_intersect(n,square_centers,first, second):
# 	f_x=square_centers[first,0]
# 	f_y=square_centers[first,1]
# 	s_x=square_centers[second,0]
# 	s_y=square_centers[second,1]

# 	if (abs(f_x-s_x)<=1.-EPS) and (abs(f_y-s_y)<=1.-EPS):
# 		return [(f_x+s_x)/2,(f_y+s_y)/2]
# 	else:
# 		return None

# def construct_P(n,square_centers):
# 	P=[]
# 	indep=np.ones(n)
# 	for i in range(n-1):
# 		for j in range(i+1,n):
# 			p=square_intersect(n,square_centers,i,j)
# 			if (p != None):
# 				P+=[p]
# 				indep[i]=indep[j]=0
# 	for i in range(n):
# 		if indep[i]>0:
# 			print(i)
# 			P+=[list(np.ravel(square_centers[i]))]
# 	return P

def is_point_in_square(p,square_centers,i):
	return (abs(p[0]-square_centers[i,0])<=0.5-EPS) and (abs(p[1]-square_centers[i,1])<=0.5-EPS)

def construct_P2(n,square_centers):
	X=np.zeros(2*n)
	Y=np.zeros(2*n)
	for i in range(n):
		X[2*i]=square_centers[i,0]-0.5
		Y[2*i]=square_centers[i,1]-0.5
		X[2*i+1]=square_centers[i,0]+0.5
		Y[2*i+1]=square_centers[i,1]+0.5
	X=np.sort(X)
	Y=np.sort(Y)
	
	Points=[]
	for i in range(2*n-1):
		for j in range(2*n-1):
			c_x=(X[i]+X[i+1])/2
			c_y=(Y[j]+Y[j+1])/2
			c=[c_x,c_y]
			inc={l for l in range(n) if is_point_in_square(c,square_centers,l)}
			if inc:
				Points.append((c,inc))
	
	lp=len(Points)

	#return Points
		
	toremove=np.zeros(lp)

	# print("{0} raw points found ...".format(lp))

	for i in range(lp-1):
		for j in range(i+1,lp):
			si=Points[i][1]
			sj=Points[j][1]
			if si.issubset(sj):
				toremove[i]=1
			elif sj.issubset(si):
				toremove[j]=1
	result=[Points[i] for i in range(lp) if toremove[i]<EPS]			
	return result



def partition_P(n,P,rotator_1):
	P_mat=np.matrix(P).transpose()
	P_rot=(rotator_1.dot(P_mat)).transpose()

	yhi=sqrt(2)/4*N
	ylow=-yhi

	lp=len(P)
	P_part=[list() for i in range(N)]

	stripe_width=sqrt(2)/2

	for i in range(lp):
		stripe=int((P_rot[i,1]-ylow)/stripe_width)
		P_part[stripe].append(i)
	return P_part

def partition_S(n,square_centers,rotator_1):
	LT=square_centers+[-0.5, 0.5]
	RT=square_centers+[ 0.5, 0.5]
	RB=square_centers+[ 0.5,-0.5]
	LB=square_centers+[-0.5,-0.5]

	yhi=sqrt(2)/4*N
	ylow=-yhi

	C_transp=(rotator_1.dot(square_centers.transpose())).transpose()
	LT_transp=(rotator_1.dot(LT.transpose())).transpose()
	RT_transp=(rotator_1.dot(RT.transpose())).transpose()
	RB_transp=(rotator_1.dot(RB.transpose())).transpose()
	LB_transp=(rotator_1.dot(LB.transpose())).transpose()

	S_part=[list() for i in range(N+1)]
	
	stripe_width=sqrt(2)/2
	for i in range(n):
		stripe_C =int((C_transp[i,1]-ylow)/stripe_width)
		# stripe_LT=int((LT_transp[i,1]-ylow)/stripe_width)
		# stripe_RT=int((RT_transp[i,1]-ylow)/stripe_width)
		# stripe_RB=int((RB_transp[i,1]-ylow)/stripe_width)
		# stripe_LB=int((LB_transp[i,1]-ylow)/stripe_width)
		# stripes={s for s in [stripe_C,stripe_LT,stripe_RT,stripe_RB,stripe_LB] if s>=0}
		# for s in stripes:
		# 	S_part[s].append(i)
		S_part[stripe_C].append(i)
	return S_part

def store_combinations(P_part):
	P_comb=[]
	P_comb_counts=[]
	for stripe in range(N):
		stripe_combinations=[]
		stripe_combinations_count=0
		Universe=P_part[stripe]
		subset_size=min(Q,len(Universe))
		for ss in range(subset_size+1):
			for C in combinations(Universe,ss):
				stripe_combinations_count+=1
				stripe_combinations.append(set(C))

		P_comb.append(stripe_combinations)
		P_comb_counts.append(stripe_combinations_count)
	return (P_comb_counts, P_comb)


def compute_T(Points,P_comb_counts,P_comb,S_part,inf):
	stripe=N-1
	# T=[]
	T=np.full((P_comb_counts[stripe-1],P_comb_counts[stripe]),inf,dtype='int32')
	squares_to_pierce=set(S_part[stripe])
	print("stripe={0}".format(stripe))
	for u_number in range(P_comb_counts[stripe-1]):
		U=P_comb[stripe-1][u_number]
		# T_U=[]
		lU=len(U)
		pierced_U=set()
		for p_number in U:
			pierced_U |= Points[p_number][1]

		for v_number in range(P_comb_counts[stripe]):
			V=P_comb[stripe][v_number]
			# print(V)
			pierced_UV=pierced_U
			# for p_number in U|V:
			for p_number in V:
				pierced_UV |= Points[p_number][1]
			if squares_to_pierce.issubset(pierced_UV):
				lUV=lU+len(V)
				T[u_number,v_number]=lUV
				# T_U.append(lUV)
			# else:
			# 	T_U.append(inf)
		# T.append(T_U)
	# print(T)
	while stripe>1:
		stripe-=1
		TT=[]
		TT=np.full((P_comb_counts[stripe-1],P_comb_counts[stripe]),inf,dtype='int32')
		print("stripe={0}".format(stripe))
		squares_to_pierce=set(S_part[stripe])
		for u_number in range(P_comb_counts[stripe-1]):
			U=P_comb[stripe-1][u_number]
			# TT_U=[]
			lU=len(U)
			pierced_U=set()
			for p_number in U:
				pierced_U |= Points[p_number][1]
			for v_number in range(P_comb_counts[stripe]):
				V=P_comb[stripe][v_number]
				TT_UV=inf
				pierced_UV=pierced_U
				for p_number in V:
					pierced_UV |= Points[p_number][1]
				for w_number in range(P_comb_counts[stripe+1]):
					W=P_comb[stripe+1][w_number]
					pierced_UVW=pierced_UV
					# for p_number in U|V|W:
					for p_number in W:
						pierced_UVW |= Points[p_number][1]
					if squares_to_pierce.issubset(pierced_UVW):
						TT_UV=min(TT_UV,lU+T[v_number][w_number])
				# TT_U.append(TT_UV)
				TT[u_number,v_number]=TT_UV
			# TT.append(TT_U)
		# print(TT)
		T=TT

	# T_final=[]
	T_final=np.full(P_comb_counts[0],inf,dtype='int32')
	squares_to_pierce=set(S_part[0])
	# print(squares_to_pierce)
	for v_number in range(P_comb_counts[0]):
		V=P_comb[0][v_number]
		T_final_V=inf
		pierced_V=set()
		for p_number in V:
			pierced_V |= Points[p_number][1]
		for w_number in range(P_comb_counts[1]):
			W=P_comb[1][w_number]
			pierced_VW=pierced_V
			# for p_number in V|W:
			for p_number in W:
				pierced_VW |= Points[p_number][1]
			if squares_to_pierce.issubset(pierced_VW):
				T_final_V=min(T_final_V,T[v_number][w_number])
		# T_final.append(T_final_V)
		T_final[v_number]=T_final_V
	return T_final


# def T(stripe,U,V,Points,P_part,S_part,inf):
# 	squares_to_pierce=set(S_part[stripe])
	
# 	# baseline
# 	if stripe == (N-1):
# 		pierced=set()
# 		for p_number in U | V:
# 			pierced |= Points[p_number][1]
# 		if squares_to_pierce.issubset(pierced):
# 			return len(U|V)
# 		else:
# 			return inf
# 	# recursion
# 	Universe=P_part[stripe+1]
# 	subset_size=min(Q,len(Universe))
# 	result=inf
	
# 	for ss in range(subset_size+1):
# 		for W in combinations(Universe,ss):
# 			pierced=set()
# 			W=set(W)
# 			for p_number in U|V|W:
# 				pierced |= Points[p_number][1]
# 			if squares_to_pierce.issubset(pierced):
# 				result=min(result, len(U)+T(stripe+1,V,W,Points,P_part,S_part,inf))
# 	return result

def Algorithm(Points,P_comb_counts, P_comb, S_part,inf):
	print("Algorithm, inf={0}".format(inf))
	T_final=compute_T(Points, P_comb_counts,P_comb,S_part,inf)
	# print(T_final)
	return min(T_final)

# def Algorithm(Points,P_part,S_part,inf):
# 	print("Algorithm, inf={0}".format(inf))
# 	Universe=P_part[0]
# 	print("Universe={0}".format(Universe))
# 	subset_size=min(Q,len(Universe))
# 	result=inf
# 	for ss in range(subset_size+1):
# 		for V in combinations(Universe,ss):
# 			print("Current V={0}".format(set(V)))
# 			T_value=T(0,set(),set(V),Points,P_part,S_part,inf)
# 			print("T-value={0}".format(T_value))
# 			result=min(result, T_value)
# 	return result


def pusk():
	start=time.time()
	if len(sys.argv)<3:
		print("No arguments, n=20 and seed=8 guessed")	
		n=20
		seed=5
	else:
		n=int(sys.argv[1])
		seed=int(sys.argv[2])

	print("{0} squares queried".format(n))
	k=-sqrt(3)
	rotator  =np.matrix([[k,1],[-1,k]])*1/sqrt(1+k**2)
	rotator_1=np.matrix([[k,-1],[1,k]])*1/sqrt(1+k**2)

	square_centers=randpoints(n,k,seed)
	# square centers (n x 2 matrix) rotated appropriately
	square_centers=(rotator.dot(square_centers.transpose())).transpose()

	print("Squares ready ...")

	Points=construct_P2(n,square_centers)
	P=[r[0] for r in Points]
	Inc=[r[1] for r in Points]
	print("{0} final points found".format(len(P)))

	#### uncomment this to employ the final Algorithm
	P_part=partition_P(n,P,rotator_1)
	S_part=partition_S(n,square_centers,rotator_1)
	# ####

	# print(P_part)
	# print("=====")
	# print(S_part)
	# print("==============")
	# print("Preprocessing complete. Proceed with the main Algorithm ...")


	#### uncomment this to employ the final Algorithm
	P_comb_counts, P_comb = store_combinations(P_part)
	####
	
	# print(P_comb_counts)	

	#### uncomment this to employ the final Algorithm
	result=Algorithm(Points, P_comb_counts, P_comb, S_part,len(Points)+1)
	####

	if result>len(Points):
		print("No piercing set")
	else:
		print("Minimum piersing set is of {0} points".format(result))

	elapsed=time.time()
	print("--- elapsed time {0:5.2f} seconds ---".format(elapsed-start))

	draw_squares(n,k,square_centers,rotator,P)


def pusk_series():
	if len(sys.argv)<2:
		print("No arguments, seed=8 guessed")	
		seed=5
	else:
		seed=int(sys.argv[1])
	for n in random.sample(range(25,300),25):
		compute_instance(n,seed)

	# draw_squares(n,k,square_centers,rotator,P)

def compute_instance(n,seed):
	start=time.time()
	# print("{0} squares queried".format(n))
	k=sqrt(3)
	rotator  =np.matrix([[k,1],[-1,k]])*1/sqrt(1+k**2)
	rotator_1=np.matrix([[k,-1],[1,k]])*1/sqrt(1+k**2)

	square_centers=randpoints(n,k,seed)
	# square centers (n x 2 matrix) rotated appropriately
	square_centers=(rotator.dot(square_centers.transpose())).transpose()

	Points=construct_P2(n,square_centers)
	P=[r[0] for r in Points]
	Inc=[r[1] for r in Points]

	# print("{0} final points found".format(len(P)))

	#### uncomment this to employ the final Algorithm	
	# P_part=partition_P(n,P,rotator_1)
	# S_part=partition_S(n,square_centers,rotator_1)

	# P_comb_counts, P_comb = store_combinations(P_part)

	# [32;115;24Mresult=Algorithm(Points,P_comb_counts, P_comb, S_part,len(Points)+1)

	# print("n={0}".format(n))
	# if result>len(Points):
	# 	print("No piercing set")
	# else:
	# 	print("Minimum piersing set is of {0} points".format(result))
	####	

	elapsed=time.time()
	# print("--- elapsed time {0:5.2f} seconds ---".format(elapsed-start))
	print("{0}, {1}".format(n,len(P)))
	sys.stdout.flush()


def com_proba():
	A=[]
	com=combinations(A,0)
	for c in com:
		print(set(c))

# com_proba()
pusk()
#pusk_series()
