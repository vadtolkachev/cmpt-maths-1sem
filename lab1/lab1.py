#!/usr/bin/python3.6

from numpy.linalg import inv
import numpy.linalg as lg
import numpy as np
import math


def multiply(matr, vec, size):
	res = np.zeros((n, 1))
	for i in range(size):
		tmp_sum = float(0)
		for j in range(size):
			tmp_sum = tmp_sum + matr[i][j]*vec[j]
		res[i] = tmp_sum

	return res


def CMInit(matr, vec, size):
	for i in range(size):
		vec[i] = 1/(i+1)
		for j in range(size):
			if(i == j):
				matr[i][j] = 1
			else:
				matr[i][j] = 1/((i+1)*(i+1) + (j+1))


def rate(x):          
	ans = 0;
	if (x.ndim > 1):
		for j in range(x.shape[1]):
			buf = 0;
			for i in range(x.shape[0]):
				buf += abs(x[i][j])
			if (buf > ans):
				ans = buf
	else:
		for j in range(x.shape[0]):
			ans += abs(x[j])
	return ans


def getDisc(A, x, f, size):
	res = multiply(A, x, size) - f
	return res


def get_Ymax(A, size):
	w, v = np.linalg.eig(A)
	yMax = np.amax(w)
	return yMax


def get_Ymin(A, size):
	w, v = np.linalg.eig(A)
	yMin = np.amin(w)
	return yMin


def swap_str(st1, st2, size):
	for i in range(size):
		st1[i], st2[i] = st2[i], st1[i]


def getMax(A, i, size):
	max_el = A[i][0]
	max_n = 0
	for j in range(size):
		if(A[i][j] > max_el):
			max_el = A[i][j]
			max_n = j
	return max_n


def getNortV(vec, size):
	res = float(0)
	for i in range(size):
		res += vec[i]*vec[i]
	
	return math.sqrt(res)
	

def getNormM(matr, size):
	norm = float(0)
	max_norm = 0
	for j in range(size):
		for i in range(size):
			norm += matr[i][j]*matr[i][j]
		if(norm > max_norm):
			max_norm = norm
		norm = 0
	return max_norm


def getU(A, size):
	u = getNormM(A, size)*getNormM(inv(A),size)
	return u


def mul(A, i, c, size):
	for j in range(size):
		A[i][j] = c*A[i][j]


def sub(A, n, m, size):
	for j in range(size):
		A[n][j] = A[n][j] - A[m][j]


def setEM(A, f, size):
	B = np.zeros((size, size+1))
	for i in range(size):
		for j in range(size):
			B[i][j] = A[i][j]
		B[i][size] = f[i]
	return B


def getLastVec(B, size):
	for i in range(size):
		f[i] = B[i][size]
	return f


def calcGauss(A, f, x, size):	
	B = setEM(A, f, size)

	for n in range(size):
		for i in range(n,size):
			mul(B, i, 1/B[i][n], size+1)
		for i in range(n+1,size):
			sub(B, i, n, size+1)
	
	for n in range(1,size):
		for j in range(size-1-n, -1, -1):
			mul(B, size-n, B[j][size-n]/B[size-n][size-n], size+1)
			sub(B, j, size-n, size+1)
			mul(B, size-n, 1/B[size-n][size-n], size+1)

	x = getLastVec(B, size)
	d = getDisc(A, x, f, size)
	return d


def mtxGauss(mtx):
	mtxabs = abs(mtx.copy())
	for j in range(mtx.shape[1] - 2):  
		maxi = np.argmax(mtxabs[j:(mtxabs.shape[0]), 0:(mtxabs.shape[1])], axis = 0)
		mtx[j], mtx[maxi[j] + j] = mtx[maxi[j] + j], mtx[j].copy()

		for i in range(mtx.shape[0]):
			if (i > j):
				mtx[i] = mtx[i] - mtx[j]*(mtx[i][j]/mtx[j][j]) 

	ans = np.zeros(mtx.shape[0])
	i = mtx.shape[0] - 1

	while (i >= 0):                         
		j = i;
		buf = mtx[i][mtx.shape[1] - 1]
		for j in range(j + 1, mtx.shape[0]):
			buf -= ans[j]*mtx[i][j]

		ans[i] = 1/mtx[i][i]*(buf)
		i -= 1
    
	return ans


def mtxZeudel(mtx, eps):

	u = np.zeros(mtx.shape[0])
    
	L = np.zeros((mtx.shape[0], mtx.shape[0]))
	D = np.zeros((mtx.shape[0], mtx.shape[0]))
	U = np.zeros((mtx.shape[0], mtx.shape[0]))
    
	for i in range(mtx.shape[0]):
		for j in range(mtx.shape[0]):
			if (i > j):
				L[i][j] = mtx[i][j]
			if (i == j):
				D[i][j] = mtx[i][j]
			if (i < j):
				U[i][j] = mtx[i][j]

	LD = lg.inv(L + D)


	k = 0
	while True: 
		uold = u
		u1 = np.dot(LD, mtx[:, mtx.shape[0]])
		u2 = -np.dot(np.dot(LD, U), u)
		u = u1 + u2
			
		if (rate(u - uold) < 2*eps): 
			break                        
		k += 1                     

	return u


n = int(12)

A = np.zeros((n, n))
f = np.zeros((n, 1))
x = np.zeros((n, 1))
d = np.zeros((n, 1))
CMInit(A, f, n)
d = calcGauss(A, f, x, n)

B = setEM(A, f, n)
f1 = np.zeros((1,n))
for i in range(n):
	f1[0][i] = 1/(i+1)

res = mtxGauss(B)
disc = (A @ res) - f1
print("Метод Гаусса\nРешение:\n", res, "\nНевязка:\n", disc, "\n")
B = setEM(A, f, n)
f1 = np.zeros((1,n))
for i in range(n):
	f1[0][i] = 1/(i+1)
e = 0.01
res = mtxZeudel(B, e)
disc = (A @ res) - f1
disc = (A @ res) - f1
print("Метод Гаусса\nРешение:\n", res, "\nНевязка:\n", disc, "\n")

Ymin = get_Ymin(A, n)
Ymax = get_Ymax(A, n)
u = getU(A, n)



print("Ymin = ", Ymin, "; Ymax = ", Ymax, ";\n")
print("u = ", u, "\n")
