import numpy as np

# input = input data. dimensions of input: q x nt, where nt = number of timesteps.
# output = response data due to input data. dimensions of output: p x nt.
# m = number of Markov parameters (impulse response timesteps) to solve for.
def okid(input,output,m=None,**options):
    input = np.atleast_2d(input)
    output = np.atleast_2d(output)

    p = output.shape[0]
    q = input.shape[0]
    nt = output.shape[1]
    assert output.shape[1] == input.shape[1]
    
    if m is None:
        m = min(300,nt)

    # adapted from Brunton
    # Form data matrix V
    V = np.zeros((q+(q+p)*m,nt))
    for i in range(nt):
        V[:q,i] = input[:q,i]
        
    for i in range(1,m+1):
        for j in range(nt-i):
            vtemp = np.concatenate((input[:,j],output[:,j]))
            V[q+(i-1)*(q+p):q+i*(q+p),i+j] = vtemp
    
    # Solve for observer Markov parameters Ybar
    Ybar = output @ np.linalg.pinv(V,rcond=10**(-3))
    
    # Isolate system Markov parameters H, and observer gain M
    D = Ybar[:,:q] # feed-through term (or D matrix) is the first term
    
    Y = np.zeros((p,q,m))
    Ybar1 = np.zeros((p,q,m))
    Ybar2 = np.zeros((p,q,m))
    
    for i in range(m):
        Ybar1[:,:,i] = Ybar[:,q+(q+p)*i : q+(q+p)*i+q]
        Ybar2[:,:,i] = Ybar[:,q+(q+p)*i+q : q+(q+p)*(i+1)]
    
    # print(Ybar2[:,:,0].shape, D.shape)
    Y[:,:,0] = Ybar1[:,:,0] + Ybar2[:,:,0] @ D
    for k in range(1,m):
        Y[:,:,k] = Ybar1[:,:,k] + Ybar2[:,:,k] @ D
        for i in range(k-1):
            Y[:,:,k] += Ybar2[:,:,i] @ Y[:,:,k-i-1]
            
    H = np.zeros((D.shape[0],D.shape[1],m+1))
    H[:,:,0] = D
    H[:,:,1:m+1] = Y[:,:,0:m]

    return H