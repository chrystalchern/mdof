import numpy as np

def construct_Hankel(data, j, start, finish):
        dim_data, nps = data.shape
        Hankel = np.empty(((finish-start+1)*dim_data, j))
        for k in range(start, finish+1): # iterate through rows
            Hankel[k*dim_data:(k+1)*dim_data,:] = data[:,start+k:start+j+k]
        return Hankel



def stacked_Hankel(inputs, outputs, j, start, finish):
        U = construct_Hankel(inputs,j,start,finish)
        Y = construct_Hankel(outputs,j,start,finish)
        # use numpy to stack U and Y vertically
        return np.vstack((U, Y)) / np.sqrt(j)# U stacked on Y divided by sqrt j


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


def extract_and_combine_blocks(matrix, row_range, col_range):
    extracted_blocks = []
    for idx in range(row_range[0], row_range[1]):
        row_blocks = []
        for jdx in range(col_range[0], col_range[1]):
            row_blocks.append(matrix[idx][jdx])
        extracted_blocks.append(row_blocks)
    new_matrix = np.block(extracted_blocks)
    return new_matrix


def choose_k(Sigma, energy_threshold=0.99):
    total_energy = np.sum(Sigma**2)
    cumulative_energy = np.cumsum(Sigma**2)
    k = np.searchsorted(cumulative_energy / total_energy, energy_threshold) + 1
    return k

def choose_k_gap(Sigma, threshold_ratio=10):
        # Calculate ratios between consecutive singular values
        ratios = Sigma[:-1] / Sigma[1:]
        # Identify the candidate k by finding the position with the maximum gap
        k_candidate = int(np.argmax(ratios)) + 1
        print("Gap ratios for each k candidate:", ratios)
        print(f"Candidate k: {k_candidate}, Maximum ratio: {ratios[k_candidate - 1]}")
        
        if ratios[k_candidate - 1] >= threshold_ratio:
            k = k_candidate
        else:
            # If no significant gap is found, consider all singular values as significant
            k = len(Sigma)
        print(f"Final selected k value: {k}")
        return k