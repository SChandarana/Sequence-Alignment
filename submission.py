import numpy as np

def dynprog(alphabet, scoring_matrix, sequence1, sequence2):
    
    def s(i,j): #function that finds out the score of a matching
        x = len(scoring_matrix)
        if i == -1: #use -1 as a parameter to match with a blank
            return scoring_matrix[x-1][alphabet.find(sequence2[j-1])]
        elif j == -1:
            return scoring_matrix[alphabet.find(sequence1[i-1])][x-1]
        else:
            return scoring_matrix[alphabet.find(sequence1[i-1])][alphabet.find(sequence2[j-1])]
    
            
    matrix = [[0 for i in range(len(sequence2)+1)] for j in range(len(sequence1)+1)] #initialise scoring matrix
    for i in range(len(sequence1)+1): #loop through it
        for j in range(len(sequence2)+1):
            if i == 0: # first column
                x = len(scoring_matrix)
                if j == 0: #top left corner
                    matrix[i][j] = scoring_matrix[x-1][x-1]
                else:
                    
                    matrix[i][j] = max(0,matrix[i][j-1] + s(-1,j))
            elif j == 0:#first row
                 matrix[i][j] = max(0,matrix[i-1][j] + s(i,-1))
            else: #everything else
                matrix[i][j] = max(0,matrix[i-1][j-1] + s(i,j), matrix[i-1][j] + s(i,-1), matrix[i][j-1] + s(-1,j)) 
    
    score = matrix[0][0]
    for i in range(len(matrix)): #acquiring indices of the highest local alignment
        if max(matrix[i]) > score:
            score = max(matrix[i])
            index_i = i
            index_j = matrix[i].index(score)

    end = False
    i = index_i
    j = index_j
    indices1 = []
    indices2 = []
    while end == False: #backtracking, end is true once i find the start of a local align
        
        if matrix[i][j] == matrix[i-1][j-1] + s(i,j):
            i -= 1
            j -= 1
            indices1.append(i)
            indices2.append(j)
        elif matrix[i][j] == matrix[i][j-1] + s(-1,j):
            j -= 1
        elif matrix[i][j] == matrix[i-1][j] + s(i,-1):
            i -= 1
        elif matrix[i][j] == 0:
            end = True
            
    indices1.reverse()
    indices2.reverse()
     
    return score, indices1, indices2

def dynproglin(alphabet, scoring_matrix, sequence1, sequence2):
    #dynprogglob is the same as dynprog but global rather than local. only called if one of the sequences is of length 1 thus linear space
    def dynprogglob(alphabet, scoring_matrix, sequence1, sequence2):
    
        def s(i,j): #function that finds out the score of a matching
            x = len(scoring_matrix)
            if i == -1: #use -1 as a parameter to match with a blank
                return scoring_matrix[x-1][alphabet.find(sequence2[j-1])]
            elif j == -1:
                return scoring_matrix[alphabet.find(sequence1[i-1])][x-1]
            else:
                return scoring_matrix[alphabet.find(sequence1[i-1])][alphabet.find(sequence2[j-1])]
        
                
        matrix = [[0 for i in range(len(sequence2)+1)] for j in range(len(sequence1)+1)]
        for i in range(len(sequence1)+1):
            for j in range(len(sequence2)+1):
                if i == 0:
                    x = len(scoring_matrix)
                    if j == 0:
                        matrix[i][j] = scoring_matrix[x-1][x-1]
                    else:
                        
                        matrix[i][j] = matrix[i][j-1] + s(-1,j)
                elif j == 0:
                    matrix[i][j] = matrix[i-1][j] + s(i,-1)
                else:
                    matrix[i][j] = max(matrix[i-1][j-1] + s(i,j), matrix[i-1][j] + s(i,-1), matrix[i][j-1] + s(-1,j)) 
        
        i = len(sequence1)
        j = len(sequence2)
        x_out = ""
        y_out = ""
        x = sequence1
        y = sequence2
        while not (i == 0 and j == 0): #returns a string of the optimal alignment
        
                
            if matrix[i][j] == matrix[i][j-1] + s(-1,j):
                j -= 1
                x_out = "-" + x_out
                y_out = y[j] + y_out
            elif matrix[i][j] == matrix[i-1][j] + s(i,-1):
                i -= 1
                y_out = "-" + y_out
                x_out = x[i] + x_out
            elif matrix[i][j] == matrix[i-1][j-1] + s(i,j):
                i -= 1
                j -= 1
                x_out = x[i] + x_out
                y_out = y[j] + y_out 

        
        return x_out, y_out
    
    def s(i,j): #function that finds out the score of a matching
        x = len(scoring_matrix)
        if i == -1: #use -1 as a parameter to match with a blank
            return scoring_matrix[x-1][alphabet.find(sequence2[j-1])]
        elif j == -1:
            return scoring_matrix[alphabet.find(sequence1[i-1])][x-1]
        else:
            return scoring_matrix[alphabet.find(sequence1[i-1])][alphabet.find(sequence2[j-1])]
    

    #NWScore and Hirschbergs are python implementations of the algorithm found here: https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
    def NWScore(x,y):
        def s(i,j): 
            temp = len(scoring_matrix)
            if i == -1: 
                return scoring_matrix[temp-1][alphabet.find(y[j-1])]
            elif j == -1:
                return scoring_matrix[alphabet.find(x[i-1])][temp-1]
            else:
                return scoring_matrix[alphabet.find(x[i-1])][alphabet.find(y[j-1])]
        
        score_matrix = [[0 for i in range(len(y)+1)],[0 for i in range(len(y)+1)]]
        
        for j in range(1,len(y)+1):
            score_matrix[0][j] = score_matrix[0][j-1] + s(-1,j)
        for i in range(1,len(x)+1):

            score_matrix[1][0] = score_matrix[0][0] + s(i,-1)
            for j in range(1,len(y)+1):
                diag = score_matrix[0][j-1] + s(i,j) 
                left = score_matrix[0][j] + s(i,-1)
                up = score_matrix[1][j-1] + s(-1,j)
                score_matrix[1][j] = max([diag,left,up])

            score_matrix[0] = list(score_matrix[1])
            
        return score_matrix[1]

    def hirschberg(x,y):
        z = ""
        w = ""
        if len(x) == 0:
            for j in range(len(y)):
                z = z + "-"
                w = w + y[j]
        elif len(y) == 0:
            for i in range(len(x)):
                z = z + x[i]
                w = w + "-"
        elif len(x) == 1 or len(y) == 1:
            (z,w) = dynprogglob(alphabet, scoring_matrix, x, y)
        else:

            xlen = len(x)
            xmid = int(len(x)/2)
            ylen = len(y)
            scoreL = NWScore(x[0:xmid],y)
            scoreR = NWScore(x[xmid:xlen][::-1],y[::-1])

            
            summed = [sum(x) for x in zip(scoreL, scoreR[::-1])]
            ymid = summed.index(max(summed))

            first = hirschberg(x[0:xmid],y[0:ymid])
            second = hirschberg(x[xmid:xlen],y[ymid:ylen])

            (z,w) = (first[0] + second[0]),(first[1] + second[1])
        
        return (z,w)
        
    #this section is used to find the start and end points of the local align in linear space
    prev_col = []
    curr_col = []
    highest_score = [-1,-1,-1]
    #loops through the two sequences and aligns them, works one column at a time so only 2 columns are needed (linear space)
    #each column stores 3 values, the score of the alignment and the i,j values of the starting point
    for i in range(len(sequence1)+1):
        prev_col = curr_col
        curr_col = []
        for j in range(len(sequence2)+1):
            if i == 0: 
                if j == 0:
                    curr_col.append([max(0,scoring_matrix[len(scoring_matrix)-1][len(scoring_matrix)-1]),i,j])
                else:
                    if 0 > curr_col[j-1][0] + s(-1,j):
                        curr_col.append([0,i,j])
                    else:
                        curr_col.append([curr_col[j-1][0] + s(-1,j),curr_col[j-1][1],curr_col[j-1][2]])
            elif j == 0:
                if 0 > prev_col[j][0] + s(i,-1):
                    curr_col.append([0,i,j])
                else:
                    curr_col.append([prev_col[j][0] + s(i,-1),prev_col[j][1],prev_col[j][2]])
            else:
                max_score = max(prev_col[j-1][0] + s(i,j), curr_col[j-1][0] + s(-1,j), prev_col[j][0] + s(i,-1))
                if 0 > max_score:
                    curr_col.append([0,i,j])
                elif max_score == prev_col[j-1][0] + s(i,j):
                    curr_col.append([max_score,prev_col[j-1][1], prev_col[j-1][2]])
                elif max_score == curr_col[j-1][0] + s(-1,j):
                    curr_col.append([max_score,curr_col[j-1][1],curr_col[j-1][2]])
                elif max_score == prev_col[j][0] + s(i,-1):
                    curr_col.append([max_score,prev_col[j][1],prev_col[j][2]])
        #print(curr_col)
        for x in range(len(curr_col)):
            if curr_col[x][0] > highest_score[0]:
                highest_score = curr_col[x]
                index_i = i
                index_j = x
    
    pointer1 = highest_score[1] #pointers used later for finding indices
    pointer2 = highest_score[2]
    local1 = sequence1[pointer1:index_i] #only contains the section of the sequence involved in the local alignment
    local2 = sequence2[pointer2:index_j]
    one,two = hirschberg(local1,local2) 
    indices1 = []
    indices2 = []
    for i in range(len(one)):
        if one[i] != "-":
            if two[i] != "-":
                indices1.append(pointer1)
                indices2.append(pointer2)
                pointer2 += 1
            pointer1 += 1
        elif two[i] != "-":
            pointer2 += 1
    



    return highest_score[0], indices1, indices2

def hueralign(alphabet, scoring_matrix, sequence1, sequence2):
    k = 2
    bandwidth = 200
    maxDiags = 10
    # 1- find all strings of length k and their indices
    # 2- find all matches of these strings and store keyed with their diagonals
    # 3- score each diagonal based on the number of seeds found within
    # 4- pick the best (maxDiags) diagonals and do bandedDP across them
    #    (starts at first seed ends at last seed) 
    # 5- return the maximum score and the alignment pattern of the best result
    def getSeeds(kval):
        seeds = {} # seed: [indices where seed occurs]
        found = False
        kval += 1
        while not found:
            kval -= 1
            for i in range(len(sequence1) + 1 - kval):
                seed = sequence1[i:i+kval]
                if seed in seeds:
                    seeds[seed].append(i)
                else:
                    seeds[seed] = [i]
            
            seedMatches = {} # diagonal value: [[i,j] of where the match occurs]
            for j in range(len(sequence2) + 1 - kval):
                toCheck = sequence2[j:j+kval]
                if toCheck in seeds:
                    found = True
                    for i in seeds[toCheck]:
                        if (i-j) in seedMatches:
                            seedMatches[i-j].append([i,j])
                        else:
                            seedMatches[i-j] = [[i,j]]
            
        return seedMatches,kval

    def createDiags(seedMatches,kval):
        diagScores = {}
        for diag in seedMatches:
            diagScores[diag] = 0
            for coords in seedMatches[diag]:
                for x in range(kval):
                    diagScores[diag] += scoring_matrix[alphabet.find(sequence1[coords[0]+x])][alphabet.find(sequence2[coords[1]+x])]
        best10scores = dict((sorted(diagScores.items(), key=lambda kv: kv[1]))[::-1][:10])
        best10diags = list(best10scores.keys())
        return best10diags   
    
    def bandedDP(diag):
        
        def s(i,j): #function that finds out the score of a matching
            x = len(scoring_matrix)
            if i == -1: #use -1 as a parameter to match with a blank
                return scoring_matrix[x-1][alphabet.find(sequence2[j-1])]
            elif j == -1:
                return scoring_matrix[alphabet.find(sequence1[i-1])][x-1]
            else:
                return scoring_matrix[alphabet.find(sequence1[i-1])][alphabet.find(sequence2[j-1])]
        
        
        k = bandwidth // 2
        if diag >= 0:
            start_i = diag - k
            start_j = 0
        else:
            start_i = 0
            start_j = -1* diag - k
            
        

        matrix = [[0 for i in range(len(sequence2)+1)] for j in range(len(sequence1)+1)] #initialise scoring matrix
        maxScore = -1
        maxScoreIndices = [-1,-1]
        
        i = start_i
        j = start_j 

        while i <= len(sequence1) + k  and j <= len(sequence2) + k :
            for x in range(j-k,j+k):
                if x >= 0 and x <= len(sequence2) and i >= 0 and i <= len(sequence1):
                    if x != 0 and i != 0:
                        diagScore = matrix[i-1][x-1] + s(i,x) 
                    else:
                        diagScore = -1
                    if x != j-k:
                        leftScore = matrix[i][x-1] + s(-1,x)
                    else:
                        leftScore = -1
                    if x != j+k:
                        upScore = matrix[i-1][x] + s(i,-1)
                    else:
                        upScore = -1
                    
                    matrix[i][x] = max(0,diagScore,leftScore,upScore)
    
                    if matrix[i][x] > maxScore:
                        maxScore = matrix[i][x]
                        maxScoreIndices = [i,x]

            i += 1
            j += 1
        
        end = False
        i = maxScoreIndices[0]
        j = maxScoreIndices[1]
        indices1 = []
        indices2 = []
        while end == False:
            if matrix[i][j] == matrix[i-1][j-1] + s(i,j):
                i -= 1
                j -= 1
                indices1.append(i)
                indices2.append(j)
            elif matrix[i][j] == matrix[i][j-1] + s(-1,j):
                j -= 1
            elif matrix[i][j] == matrix[i-1][j] + s(i,-1):
                i -= 1
            elif matrix[i][j] == 0:
                end = True
                
        indices1.reverse()
        indices2.reverse()
        return maxScore, indices1, indices2
            
                    
    maxScore = -1
    indices1 = []
    indices2 = []
    
    seedMatches,kval = getSeeds(k)
    
    best = createDiags(seedMatches,kval)
    for diag in best:
        tempMax, tempInd1, tempInd2 = bandedDP(diag)
        if tempMax > maxScore:
            maxScore = tempMax
            indices1 = tempInd1
            indices2 = tempInd2
    return maxScore, indices1, indices2
