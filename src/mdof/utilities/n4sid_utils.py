import numpy as np



u = np.arange(1,121).reshape((40,3)).T
y = np.arange(1,41).reshape((40,1)).T


def construct_hankel(data, j, start, finish):
    dim_data, nps = data.shape
    Hankel = np.empty(((finish-start+1)*dim_data, j))
    for k in range(start, finish+1):
        Hankel[k*dim_data:(k+1)*dim_data,:] = data[:, start+k:start+j+k]
    return Hankel

def stacked_hankel(u, y, j, start, finish):
    U = construct_hankel(u, j, start, finish)
    Y = construct_hankel(y, j, start, finish)
    return np.vstack((U, Y)) / np.sqrt(j)

def slice_matrix(matrix, row_sizes, col_sizes):
        blocks = []
        row_start = 0
        for row_size in row_sizes:
            col_start = 0
            row_end = row_start + row_size
            row_blocks = []
            for col_size in col_sizes:
                col_end = col_start + col_size
                block = matrix[row_start:row_end, col_start:col_end]
                row_blocks.append(block)
                col_start += col_size
            blocks.append(row_blocks)
            row_start += row_size
        return blocks


def partition_R_matrix(RT, QT, i, j, m, l):
    RT_row_sizes = [3*i, 3, 3*(i-1), i, 1, i-1] 
    RT_col_sizes = [3*i, 3, 3*(i-1), i, 1, i-1] 
    QT_row_sizes = [3*i, 3, 3*(i-1), i, 1, i-1]  
    QT_col_sizes = [j]  

    RT_blocks = slice_matrix(RT, RT_row_sizes, RT_col_sizes)
    QT_blocks = slice_matrix(QT, QT_row_sizes, QT_col_sizes)
    
    return RT_blocks, QT_blocks

def compute_projection_matrices(RT_blocks, R_blocks, i, l, m):
    R5614 = [[RT_blocks[i][j] for j in range(4)] for i in range(4, 6)]
    new_matrix_R5614 = np.block(R5614)

    RT1414 = [[R_blocks[i][j] for j in range(4)] for i in range(4)]
    new_matrix_RT1414 = np.block(RT1414)

    return new_matrix_R5614, new_matrix_RT1414

def compute_Li_matrices(new_matrix_R5614, new_matrix_RT1414, i, l, m):
    result_matrix = np.dot(new_matrix_R5614, new_matrix_RT1414)

    l_i, m_i = l * i, m * i
    Li1 = result_matrix[:l_i, :m_i]
    Li2 = result_matrix[:l_i, m_i:2*m_i]
    Li3 = result_matrix[:l_i, 2*m_i:]

    return Li1, Li2, Li3

def compute_gamma_and_svd(Li1, Li3, i, l, m):
    Zero_Matrix_1 = np.zeros((l*i, m*i))
    Combined_Matrix = np.hstack((Li1, Zero_Matrix_1, Li3))
    
    U, Sigma, VT = np.linalg.svd(Combined_Matrix, full_matrices=False)
    
    k = 2  
    U1 = U[:, :k]
    Sigma1 = np.diag(Sigma[:k])
    
    Sigma1_half = np.sqrt(Sigma1)
    Gamma_i = np.dot(U1, Sigma1_half)

    return Gamma_i, U1, Sigma1

def compute_state_space_matrices(U1, Sigma1, RT_blocks, R_blocks, i, l, m):
    Sigma1_neg_half = np.diag(1 / np.sqrt(np.diag(Sigma1)))
    
    R1514 = [[RT_blocks[i][j] for j in range(4)] for i in range(5)]
    new_matrix_R1514 = np.block(R1514)
    
    R5514 = [[RT_blocks[i][j] for j in range(4)] for i in range(4, 5)]
    new_matrix_R5514 = np.block(R5514)

    R2214 = [[RT_blocks[i][j] for j in range(4)] for i in range(1, 2)]
    new_matrix_R2214 = np.block(R2214)
    
    final_matrix_3 = np.vstack((np.dot(Sigma1_neg_half, U1.T), new_matrix_R5514))
    final_matrix_4 = np.vstack((np.dot(Sigma1_neg_half, U1.T), new_matrix_R2214))

    AAA_pinv = np.linalg.pinv(final_matrix_4)
    script_L = np.dot(final_matrix_3, AAA_pinv)

    A, B = script_L[:i, :i], script_L[:i, i:i+m]
    C, D = script_L[i:i+l, :i], script_L[i:i+l, i:i+m]

    
    k=2
    j=26
    rho_2 = final_matrix_3 - np.dot(script_L, final_matrix_4)
    print(rho_2)

    rho1_2 = rho_2[:k, :]
    rho2_2 = rho_2[k:, :]
    print("rho1_2 (first k rows of rho_2):")
    print(rho1_2)
    print("\nrho2_2 (rows from k+1 to the end of rho_2):")
    print(rho2_2)

    rho1_2_transposed = rho1_2.T
    rho2_2_transposed = rho2_2.T

    Qs = (1/j) * np.dot(rho1_2, rho1_2_transposed)
    Ss = (1/j) * np.dot(rho1_2, rho2_2_transposed)
    Rs = (1/j) * np.dot(rho2_2, rho2_2_transposed)
    print("Qs =")
    print(Qs)
    print("\nSs =")
    print(Ss)
    print("\nRs =")
    print(Rs)

    system = (A, B, C, D)
    from mdof.simulate import simulate
    output_prediction = simulate(system, u)
    import matplotlib.pyplot as plt
    plt.plot(y.T)
    print(y)
    plt.show()
    print(output_prediction.T)
    plt.plot(output_prediction.T)
    plt.show()

    my_sigma = np.array([5,4,3,2,1,1e-10,1e-11])

    return A, B, C, D