#!/usr/bin/env python3

"""
Author: 

Description: this is a script to ...
"""
#import statements here

# functions between here and __main__
blosum = """
# http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
   A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
   R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
   N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
   D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
   C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
   Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
   E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
   H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
   I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
   L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
   K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
   M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
   F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
   P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
   S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
   T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
   W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
   Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
   V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
   B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
   Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
   X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
   * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""

def blosum62():
    """Return order and similarity scores from BLOSUM62 matrix

    order: dict of {res: idx_in_matrix}
    blosum_matrix: list of lists with similarity scores
    """
    order = {}
    blosum_matrix = []
    for line in blosum.split('\n'):
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) == 24:
            for idx, sym in enumerate(parts):
                order[sym] = idx
        else:
            # list around the map construction for python3 compatibility
            blosum_matrix.append(list(map(int,parts[1:])))
    return order, blosum_matrix

BLOSUM62_ORDER, BLOSUM62_MATRIX = blosum62()

def score(res1, res2, BLOSUM62_ORDER, BLOSUM62_MATRIX):
    """Return similarity score from BLOSUM62 matrix for two residues
    
    res1: string, amino acid
    res2: string, amino acid
    """
    lookup1 = BLOSUM62_ORDER[res1]
    lookup2 = BLOSUM62_ORDER[res2]
    return BLOSUM62_MATRIX[lookup1][lookup2]

# write your own functions below here

def decide_score(input_matrix,x,y,seq1,seq2,endgap_bool = False,gap_penalty = -8,end_gap_penalty = -8): #the gap penalty function
    """Return similarity score and optimal route direction
    
    input_matrix: 2d array, score matrix
    x,y: int, coordinates
    seq1,seq2: string, the strings whose residues are scored
    endgap_bool: boolean, wether an end gap is applied
    gap_penalty: int, the gap penalty
    end_gap_penalty: int, the end gap penalty
    """

    if endgap_bool:
        return end_gap_penalty
    
    ho = input_matrix[y][x-1] + gap_penalty #horizontal

    ve =  input_matrix[y-1][x] + gap_penalty #vertical
        
    di =  input_matrix[y-1][x-1] + score(seq1[y-1],seq2[x-1]) #diagonal
    
    if max(di,ho,ve) == di:
        return di,'di'
    elif max(di,ho,ve) == ho:
        return ho,'ho'
    else:
        return ve,'ve'
        
    
def NW_matrix(seq1, seq2,endgap_bool = False,penalty = -8, endpenalty = -8):
    """Returns 2 matrices, the score matrix and the direction matrix
    
    seq1: string, first sequence to be compared. Y axis
    seq2: string, second sequence to be compared. X axis
    endgap_bool: boolean, whether there is an endgap
    penalty: int, penalty amount
    endpenalty: int, endgap penalty amount
    """

    output_matrix = [[0 for x in range(len(seq2)+1)] for y in range(len(seq1)+1)]
    direction_matrix = [[0 for x in range(len(seq2)+1)] for y in range(len(seq1)+1)]
    output_matrix[0][0] = 0 

    if endgap_bool:#initialize endgaps if applicable
        for i in range(1,len(output_matrix)):
            output_matrix[i][0] = output_matrix[i-1][0] + decide_score(output_matrix,0,i,seq1,seq2,True,penalty,endpenalty)

        for j in range(1,len(output_matrix[0])):
            output_matrix[0][j] = output_matrix[0][j-1]+ decide_score(output_matrix,j,0,seq1,seq2,True,penalty,endpenalty)

    for y in range(1,len(output_matrix)): #y iterator
        
        for x in range(1,len(output_matrix[0])): #x iterator
            #assign the matrix element the correct score for its current position
            #print(str(y) + str(x))
            output_matrix[y][x],direction_matrix[y][x] = decide_score(output_matrix,x,y,seq1,seq2,gap_penalty = penalty,end_gap_penalty = endpenalty)
            
                
        #print('y: ' + str(len(output_matrix)) + ' ' + 'x: ' +str(len(output_matrix[0])))
    return output_matrix,direction_matrix

        
def traceback(matrix):
    """Return a path of the optimal alignment direction through the score matrix, starting at the end and going to the beginning
    matrix: 2d array, a matrix of direction arguments saved as strings
    """
    path = [] #saves the path,0 for diagonal, 1 for horizontal, 2 for vertical
    x = len(matrix[0]) -1
   
    y = len(matrix) -1
    
    while x >= 0 and y >= 0:
        if matrix[y][x] == 'di':
            if path:
                path.append('di')
            else:
                path = ['di']
            x -= 1
            y -= 1
        elif matrix[y][x] == 'ho':
            if path:
                path.append('ho')
            else:
                path = ['ho']
            x-= 1
        else:
            if path:
                path.append('ve')
            else:
                path = ['ve']
            y-= 1
        
        
    return path

def visualize_alignment(seq1,seq2,path):
    """Returns nothing
    visualizes the matrix, through interpreting the path, prints results
    seq1: string, first sequence to be compared. Y axis
    seq2: string, second sequence to be compared. X axis
    path: list, the directional optimal path through the matrix
    """
    mPath = path[::-1]
    string_seq1 = ''''''
    iterator_seq1 = 0
    hasmoved_vertical = False
    string_seq2 = ''''''
    iterator_seq2 = 0
    hasmoved_horizontal = False
    
    
    #print(len(path))
    #path code reminder: 0 for diagonal, 1 for horizontal, 2 for vertical
    for direction in mPath:
        if direction == 'di':#diagonal, both strings gain a piece of the sequence
            
            string_seq1 += seq1[iterator_seq1]
            iterator_seq1 +=1

            string_seq2 += seq2[iterator_seq2]
            iterator_seq2 +=1
        
            hasmoved_vertical = True
            hasmoved_horizontal = True
            #set alignment
            #string_alignment += '|'
        
                

        if direction == 'ho':#horizontal, gap in Y, x sequence gains a piece
            if hasmoved_horizontal:
                string_seq2 += seq2[iterator_seq2]
                iterator_seq2 +=1
                #string_alignment += ' '

            string_seq1 += '-'
            hasmoved_horizontal = True

        if direction == 've':#vertical, gap in X, y sequence gains a piece
            if hasmoved_vertical:
                string_seq1 += seq1[iterator_seq1]
                iterator_seq1 +=1
                string_seq2 += '-'
                #string_alignment += ' '
            hasmoved_vertical = True
        #set alignment better
        string_alignment = ''''''
        for i in range(len(string_seq1)):
            if string_seq1[i] == string_seq2[i]:
                string_alignment += '|'
            else:
                string_alignment += ' '

            
    
    print('\nAlignment:\n\n' + string_seq1 + '\n' + string_alignment + '\n' + string_seq2)
    #print(str(len(string_seq1)) + '\n' + str(len(string_seq2)) + '\n' + str(len(string_alignment)))
    permatch = percentage_match(string_alignment)
    print("{0:.2f}%".format(percentage_match(string_alignment)) + ' Match')

def percentage_match(alignment):
    """Returns output percentage of alignment between sequences
    
    alignment: list, the alignment given in | and spaces
    """
    matches = 0
    for c in alignment:
        if c == '|':
            matches += 1

    output_percentage = float(matches / len(alignment) * 100 )
    return output_percentage
    

def align_sequences(seq1,seq2,endgap_bool = False,penalty = -8, endpenalty = -8):
    """Returns nothing
    
    aligns two sequences and prints the relevant information to the console

    seq1: string, first sequence to be compared. Y axis\n
    seq2: string, second sequence to be compared. X axis\n
    endgap_bool: boolean, whether there is an endgap\n
    penalty: int, penalty amount\n
    endpenalty: int, endgap penalty amount
    """
    print('aligning sequences with gap: ' + str(endgap_bool) + '\nPenalty = ' + str(penalty))
    print('sequence 1 = ' + seq1 + '\nsequence 2 = ' + seq2 )
    score_matrix,direction_matrix = NW_matrix(seq1,seq2,endgap_bool,penalty,endpenalty)
    print('\nMatrix:\n')
    for line in score_matrix:
        print(str(line))
    traceback_path = traceback(direction_matrix)
    
    visualize_alignment(seq1,seq2,traceback_path)


