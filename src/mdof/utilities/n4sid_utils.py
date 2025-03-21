import numpy as np

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

'''
def partition_R_matrix(RT, QT, i, j, m, l):
    RT_row_sizes = [m*i, m, m*(i-1), l*i, l, l*(i-1)] 
    RT_col_sizes = [m*i, m, m*(i-1), l*i, l, l*(i-1)] 
    QT_row_sizes = [m*i, m, m*(i-1), l*i, l, l*(i-1)]  
    QT_col_sizes = [j]  

    RT_blocks = slice_matrix(RT, RT_row_sizes, RT_col_sizes)
    QT_blocks = slice_matrix(QT, QT_row_sizes, QT_col_sizes)
    
    return RT_blocks, QT_blocks

def partition_R1_matrix(R, i, j, m, l):
    R_row_sizes = [m*i, m, m*(i-1), l*i, l, l*(i-1)] 
    R_col_sizes = [m*i, m, m*(i-1), l*i, l, l*(i-1)] 
    R_blocks = slice_matrix(R, R_row_sizes, R_col_sizes)
    
    return R_blocks
'''

def partition_R_matrices(R, RT, i, j, m, l):
    row_sizes = [m*i, m, m*(i-1), l*i, l, l*(i-1)] 
    col_sizes = [m*i, m, m*(i-1), l*i, l, l*(i-1)]
    RT_blocks = slice_matrix(RT, row_sizes, col_sizes)
    R_blocks = slice_matrix(R, row_sizes, col_sizes)

    return R_blocks, RT_blocks

def compute_projection_matrices(RT_blocks,R_blocks):
    R5614 = [[RT_blocks[i_idx][j_idx] for j_idx in range(4)] for i_idx in range(4, 6)]
    new_matrix_R5614 = np.block(R5614)

    # TODO: recompute this part with inverse of R
    RT1414 = [[R_blocks[i_idx][j_idx] for j_idx in range(4)] for i_idx in range(4)]
    new_matrix_RT1414 = np.block(RT1414)

    R6615 = [[RT_blocks[i_idx][j_idx] for j_idx in range(5)] for i_idx in range(5,6)]
    new_matrix_R6615 = np.block(R6615)

    RT1515 = [[R_blocks[i_idx][j_idx] for j_idx in range(5)] for i_idx in range(5)]
    new_matrix_RT1515 = np.block(RT1515)

    return new_matrix_R5614, new_matrix_RT1414, new_matrix_R6615, new_matrix_RT1515

def compute_Li_matrices(new_matrix_R5614, new_matrix_RT1414, i, l, m):
    result_matrix = np.dot(new_matrix_R5614, new_matrix_RT1414)

    l_i, m_i = l * i, m * i
    Li1 = result_matrix[:l_i, :m_i]
    Li2 = result_matrix[:l_i, m_i:2*m_i]
    Li3 = result_matrix[:l_i, 2*m_i:]

    return Li1, Li2, Li3

def compute_Li_1_matrices(new_matrix_R6615, new_matrix_RT1515, i, l, m):
    result_matrix_1 = np.dot(new_matrix_R6615, new_matrix_RT1515)

    m_i = m * i
    Li_11 = result_matrix_1[:l*(i-1), :m*(i+1)] 
    Li_12 = result_matrix_1[:l*(i-1), m*(i+1):2*m_i]  
    Li_13 = result_matrix_1[:l*(i-1), 2*m_i:]

    return Li_11, Li_12, Li_13

def compute_gamma_and_svd(Li1, Li3, i, l, m, RT_blocks, threshold=1e-3):
    Zero_Matrix_1 = np.zeros((l*i, m*i))
    Combined_Matrix = np.hstack((Li1, Zero_Matrix_1, Li3))
    R1414 = [[RT_blocks[i_idx][j_idx] for j_idx in range(4)] for i_idx in range(4)]
    new_matrix_R1414 = np.block(R1414)
    result_matrix_2 = np.dot(Combined_Matrix, new_matrix_R1414)
    U, Sigma, VT = np.linalg.svd(result_matrix_2, full_matrices=False)
    
    
    k = int(np.sum(Sigma > threshold))
    
    U1 = U[:, :k]
    Sigma1 = np.diag(Sigma[:k]) 
    Sigma1_half = np.sqrt(Sigma1)
    Gamma_i = np.dot(U1, Sigma1_half)

    return Gamma_i, U1, Sigma1, k

def compute_state_space_matrices(U1, Sigma1, RT_blocks, i, l, m, Li1, Li3, Li_11, Li_13, k, j, verbose=False, stochastic_terms=True):
    Sigma1_neg_half = np.diag(1 / np.sqrt(np.diag(Sigma1)))
    
    U1_truncated = U1[:-l, :]
    U1_truncated_pinv = np.linalg.pinv(U1_truncated)

    R1514 = [[RT_blocks[i_idx][j_idx] for j_idx in range(4)] for i_idx in range(5)]
    new_matrix_R1514 = np.block(R1514)
    
    R5514 = [[RT_blocks[i_idx][j_idx] for j_idx in range(4)] for i_idx in range(4, 5)]
    new_matrix_R5514 = np.block(R5514)

    R2214 = [[RT_blocks[i_idx][j_idx] for j_idx in range(4)] for i_idx in range(1, 2)]
    new_matrix_R2214 = np.block(R2214)

    R1414 = [[RT_blocks[i_idx][j_idx] for j_idx in range(4)] for i_idx in range(4)]
    new_matrix_R1414 = np.block(R1414)
    
    Zero_Matrix_1 = np.zeros((l*i, m*i))
    Combined_Matrix = np.hstack((Li1, Zero_Matrix_1, Li3))
    
    Zero_Matrix_2 = np.zeros((l*(i-1), m*(i-1)))
    Combined_Matrix_2 = np.hstack((Li_11, Zero_Matrix_2, Li_13))
    
    final_matrix_3 = np.vstack((Sigma1_neg_half @ U1_truncated_pinv @ Combined_Matrix_2 @ new_matrix_R1514, new_matrix_R5514))
    final_matrix_4 = np.vstack((Sigma1_neg_half @ U1.T @ Combined_Matrix @ new_matrix_R1414, new_matrix_R2214))
    

    AAA_pinv = np.linalg.pinv(final_matrix_4)
    script_L = np.dot(final_matrix_3, AAA_pinv)

    A, B = script_L[:k, :k], script_L[:k, k:k+m]
    C, D = script_L[k:k+l, :k], script_L[k:k+l, k:k+m]
    
    Qs, Ss, Rs = None, None, None

    if stochastic_terms:
        rho_2 = final_matrix_3 - np.dot(script_L, final_matrix_4)

        rho1_2 = rho_2[:k, :]
        rho2_2 = rho_2[k:, :]

        rho1_2_transposed = rho1_2.T
        rho2_2_transposed = rho2_2.T

        Qs = (1/j) * np.dot(rho1_2, rho1_2_transposed)
        Ss = (1/j) * np.dot(rho1_2, rho2_2_transposed)
        Rs = (1/j) * np.dot(rho2_2, rho2_2_transposed)

    return A,B,C,D,Qs,Ss,Rs
    