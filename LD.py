phi(r):
	return 4*epsilon * ((sigma/r)**12 - (sigma/r)**6)

phi_dash(r):
	return epsilon*((24*sigma**6/r**7) - (48*sigma**12/r**13))
	
phi_ddash(r):
	return epsilon*((624*sigma**12/r**14) - (168*sigma**6/r**8))



findX(j,k,jdash,kdash,alpha,beta):
	r_alpha = Pos[j,k,alpha]- Pos[jdash,kdash,alpha]
	r_beta = Pos[j,k,beta]- Pos[jdash,kdash,beta] 
	r = ((Pos[j,k,0]- Pos[jdash,kdash,0])**2 + (Pos[j,k,1]- Pos[jdash,kdash,1])**2 + (Pos[j,k,2]- Pos[jdash,kdash,2])**2)**0.5

	X = (r_alpha * r_beta /r**2)*(phi_ddash(r)-(1/r*phi_dash(r)))
	
	if alpha != beta:
		X = X + phi_dash(r)/r
		
	return (r)

## Generate the 3x3 force constant matrix ##

forceConstant(j,k,jdash, kdash):
	
	phi = np.zeros((3,3))	
	for alpha in range(3):
		for beta in range(3):
			if j == jdash and k == kdash:
				X = 0
				for kddash in range(k):
					for jddash in range(n):
						X = X + findX(j,k,jddash,kddash,alpha,beta)
				phi[alpha,beta] = X

			else:
				X = findX(j,k,jdash,kdash,alpha,beta)
				phi[alpha,beta] = -X
									

