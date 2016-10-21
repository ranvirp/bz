import numpy
import thinning
import itertools as it

def guo_hall_thinning(src):
	def iteration(src, iter):
		marker = numpy.ones(src.shape, numpy.uint8)
		h,w = src.shape
		changed = 0
		for j,i in numpy.transpose(numpy.nonzero(src)):
			if i==0 or i==w-1: continue
			if j==0 or j==h-1: continue
			assert src.item(j,i)!=0
			p2 = src.item((j,   i-1))
			p3 = src.item((j+1, i-1))
			p4 = src.item((j+1, i))
			p5 = src.item((j+1, i+1))
			p6 = src.item((j,   i+1))
			p7 = src.item((j-1, i+1))
			p8 = src.item((j-1, i))
			p9 = src.item((j-1, i-1))
			C = ((~p2 & (p3 | p4)) + 
				(~p4 & (p5 | p6)) + 
				(~p6 & (p7 | p8)) +
				(~p8 & (p9 | p2)))
			N1 = (p9 | p2) + (p3 | p4) + (p5 | p6) + (p7 | p8)
			N2 = (p2 | p3) + (p4 | p5) + (p6 | p7) + (p8 | p9)
			N = min(N1, N2)
			if iter==0:
				m = (p8 & (p6 | p7 | ~p9))
			else:
				m = (p4 & (p2 | p3 | ~p5))
			# print i,j
			# print p2, p3,p4,p5,p6,p7,p8,p9
			# print "C",C, "N",N, "m", m
			if C==1 and 2<= N <=3 and m==0:
				marker.itemset((j,i),0)
				changed += 1
		return src & marker, changed

	dst = src.copy()
	i=0;
	while True:
		i+=1
		# now = time.clock()
		dst, changed  = iteration(dst, 0)
		dst, changed2 = iteration(dst, 1)
		# dst, changed  = c_iteration(dst, 0)
		# dst, changed2 = c_iteration(dst, 1)

		# print time.clock() - now
		d = changed + changed2
		if d == 0:
			break
			
	return dst 

def test():
	i=0
	for w in range(1,5):
		for h in range(1,5):
			arr = numpy.empty(w*h, dtype=numpy.ubyte)
			for o in it.product([0,1],repeat=w*h):
				arr.flat = o
				arr = numpy.reshape(arr,[h,w])
				print i
				# print arr
				i+=1
				# print arr
				# arr[1][1]=0
				brr = thinning.guo_hall_thinning(arr.copy())
				arr2 = guo_hall_thinning(arr)
				if not all((arr2 == brr).flat):
					print "original:",arr
					print "c",brr
					print "python", arr2
					return 
			# print brr

test()