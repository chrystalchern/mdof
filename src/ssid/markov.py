import numpy as np

# y = output data: forced response. dimensions of y: p x nt, where nt = number of timesteps
# u = input data: input. dimensions of u: q x nt.
# m = number of Markov parameters (impulse response timesteps) to solve for
def okid(y,u,m):
    y = np.atleast_2d(y)
    u = np.atleast_2d(u)
    
    # adapted from Brunton
    p = y.shape[0]
    q = u.shape[0]
    nt = y.shape[1]
    assert y.shape[1] == u.shape[1]
    
    # Form data matrices y and V
    V = np.zeros((q+(q+p)*m,nt))
    for i in range(nt):
        V[:q,i] = u[:q,i]
        
    for i in range(1,m+1):
        for j in range(nt-i):
            vtemp = np.concatenate((u[:,j],y[:,j]))
            V[q+(i-1)*(q+p):q+i*(q+p),i+j] = vtemp
    
    # Solve for observer Markov parameters Ybar
    Ybar = y @ np.linalg.pinv(V,rcond=10**(-3))
    
    # Isolate system Markov parameters H, and observer gain M
    D = Ybar[:,:q] # feed-through term (or D matrix) is the first term
    
    Y = np.zeros((p,q,m))
    Ybar1 = np.zeros((p,q,m))
    Ybar2 = np.zeros((p,q,m))
    
    for i in range(m):
        Ybar1[:,:,i] = Ybar[:,q+(q+p)*i : q+(q+p)*i+q]
        Ybar2[:,:,i] = Ybar[:,q+(q+p)*i+q : q+(q+p)*(i+1)]
    
    Y[:,:,0] = Ybar1[:,:,0] + Ybar2[:,:,0] @ D
    for k in range(1,m):
        Y[:,:,k] = Ybar1[:,:,k] + Ybar2[:,:,k] @ D
        for i in range(k-1):
            Y[:,:,k] += Ybar2[:,:,i] @ Y[:,:,k-i-1]
            
    H = np.zeros((D.shape[0],D.shape[1],m+1))
    H[:,:,0] = D
    
    for k in range(1,m+1):
        H[:,:,k] = Y[:,:,k-1]
        
    return H


def OKID(y,u,r):
    # inputs:  y (sampled output), u (sampled input), r (effective system order)
    # outputs: H (Markov parameters), M (Observer gain)
    
    PP = y.shape[0] # number of outputs
    MM = y.shape[1] # number of output samples
    QQ = u.shape[0] # number of inputs
    lu = u.shape[1] # number of input samples
    
    assert MM == lu
    print(MM)
    print(lu)
    
    LL = r*5
    
    # Form data matrices y and V
    V = np.zeros((QQ+(QQ+PP)*LL,MM))
    for i in range(MM):
        V[:QQ,i] = u[:QQ,i]
        
    for i in range(1,LL+1):
        for j in range(MM-i):
            vtemp = np.concatenate((u[:,j],y[:,j]))
            V[QQ+(i-1)*(QQ+PP):QQ+i*(QQ+PP),i+j] = vtemp
    
    # Solve for observer Markov parameters Ybar
    Ybar = y @ np.linalg.pinv(V,rcond=10**(-3))
    
    # Isolate system Markov parameters H, and observer gain M
    D = Ybar[:,:QQ] # feed-through term (or D matrix) is the first term
    
    Y = np.zeros((PP,QQ,LL))
    Ybar1 = np.zeros((PP,QQ,LL))
    Ybar2 = np.zeros((PP,QQ,LL))
    
    for i in range(LL):
        Ybar1[:,:,i] = Ybar[:,QQ+(QQ+PP)*i : QQ+(QQ+PP)*i+QQ]
        Ybar2[:,:,i] = Ybar[:,QQ+(QQ+PP)*i+QQ : QQ+(QQ+PP)*(i+1)]
    
    Y[:,:,0] = Ybar1[:,:,0] + Ybar2[:,:,0] @ D
    for k in range(1,LL):
        Y[:,:,k] = Ybar1[:,:,k] + Ybar2[:,:,k] @ D
        for i in range(k-1):
            Y[:,:,k] += Ybar2[:,:,i] @ Y[:,:,k-i-1]
            
    H = np.zeros((D.shape[0],D.shape[1],LL+1))
    H[:,:,0] = D
    
    for k in range(1,LL+1):
        H[:,:,k] = Y[:,:,k-1]
        
    return H