import pandas as pandas
#test

def init(rows, cols, gap_penalty=10):
    r = []
    c = []
    for j in range(0,cols+1):
        c.append(j*-1*gap_penalty)
    r.append(c)
    for i in range(1, rows+1):
        c = [i*-1*gap_penalty]
        for j in range(1,cols+1):
            c.append(0.0)
        r.append(c)
    return r

def get_new_score(up, left, middle,
                  matched, gap_penalty, match, mismatch):
    return max([middle + (match if matched == 1 else mismatch), up - gap_penalty, left - gap_penalty])


def get_alignment_u(top_seq, bottom_seq, sm,
                  gap_penalty, match, mismatch):
    aligned1 = ''
    aligned2 = ''
    i = len(top_seq)
    j = len(bottom_seq)
    while i > 0 or j > 0:
        score = sm[i][j]
        score_diag = sm[i-1][j-1]
        score_left = sm[i][j-1]
        score_up = sm[i-1][j]
        if i > 0 and j > 0 and score_diag ==  max(score_diag,score_left,score_up):
            aligned1 = top_seq[i-1] + aligned1
            aligned2 = bottom_seq[j-1] + aligned2
            i -= 1
            j -= 1
        else:
            if j > 0 and score_left ==  max(score_diag,score_left,score_up):
                aligned1 = '-' + aligned1
                aligned2 = bottom_seq[j - 1] + aligned2
                j -= 1
            else:
                #if score == score_left - gap_penalty:
                aligned1 = top_seq[i - 1] + aligned1
                aligned2 = '-' + aligned2
                i -= 1
    return "{}\n{}".format(aligned1, aligned2)

def align_u(top_seq, bottom_seq, gap_penalty=10,
          match=2, mismatch=-1):
    top_len = len(top_seq)
    bottom_len = len(bottom_seq)
    mrx = init(top_len, bottom_len,gap_penalty)
    for i in range(1,top_len+1):
        for j in range(1,bottom_len+1):
            mrx[i][j] = get_new_score(mrx[i-1][j],
                                      mrx[i][j-1],
                                      mrx[i-1][j-1],
                                      1 if top_seq[i-1] == bottom_seq[j-1] else 0,
                                      gap_penalty,
                                      match,
                                      mismatch
                                      )
    #for row in mrx:
    #    print(' '.join(map(str, row)))
    #import pandas
    print(pandas.DataFrame(mrx))
    aligns = get_alignment_u(top_seq, bottom_seq, mrx, gap_penalty, match, mismatch)
    return aligns, mrx[top_len][bottom_len]

def get_alignment(top_seq, bottom_seq, sm,
                  gap_penalty, match, mismatch):
    aligned1 = ''
    aligned2 = ''
    i = len(top_seq)
    j = len(bottom_seq)
    while i > 0 and j > 0:
        score = sm[i][j]
        score_diag = sm[i-1][j-1]
        score_left = sm[i][j-1]
        score_up = sm[i-1][j]
        if score == score_diag + (match if top_seq[i-1] == bottom_seq[j-1] else mismatch):
            aligned1 = top_seq[i-1] + aligned1
            aligned2 = bottom_seq[j-1] + aligned2
            i -= 1
            j -= 1
        else:
            if score == score_left - gap_penalty:
                aligned1 = '-' + aligned1
                aligned2 = bottom_seq[j - 1] + aligned2
                j -= 1
            else:
                if score ==  - gap_penalty:
                    aligned1 = top_seq[i - 1] + aligned1
                    aligned2 = '-' + aligned2
                    i -= 1

    while i > 0 :
        aligned1 = top_seq[i-1] + aligned1
        aligned2 = '-' + aligned2
        i -= 1
    while j > 0 :
        aligned1 = '+' - aligned1
        aligned2 = bottom_seq[get_alignment_u()-1] + aligned2
        j -= 1

    return "{}\n{}".format(aligned1, aligned2)

def align(top_seq, bottom_seq, gap_penalty=10,
          match=2, mismatch=-1):
    top_len = len(top_seq)
    bottom_len = len(bottom_seq)
    mrx = init(top_len, bottom_len,gap_penalty)
    for i in range(1,top_len+1):
        for j in range(1,bottom_len+1):
            mrx[i][j] = get_new_score(mrx[i-1][j],
                                      mrx[i][j-1],
                                      mrx[i-1][j-1],
                                      1 if top_seq[i-1] == bottom_seq[j-1] else 0,
                                      gap_penalty,
                                      match,
                                      mismatch
                                      )
    #for row in mrx:
    #    print(' '.join(map(str, row)))
    #import pandas
    #print(pandas.DataFrame(mrx))
    aligns = get_alignment(top_seq, bottom_seq, mrx, gap_penalty, match, mismatch)
    return aligns, mrx[top_len][bottom_len]




#a,s = align('GCATGCU', 'GATTACA', 1, 1, -1)
#print(a, '\n', s, '\n')

#a,s = align('ACTTCTACCCCAGC', 'AGTGTCGGCT', 1.399, 2.2168, -4.4499)
#print(a, '\n', s, '\n')

##a,s = align('ACTTCTACCCCAGC', 'AGTGTCGGCT', 1.399, 2.2168, -4.4499)
##print(a, '\n', s, '\n')

a,s = align('AAGATTTTAGATT', 'GCTTCCTTTTCGCGAGCGT', 1.399, 2.2168, -4.4499)
print(a, '\n', s, '\n')

#a,s = align('AGTGTCGGCT', 'ACTTCTACCCCAGC', 1.399, 2.2168, -4.4499)
#print(a, '\n', s, '\n')

#a,s = align_u('AGTGTCGGCT', 'ACTTCTACCCCAGC', 1.399, 2.2168, -4.4499)
#print(a, '\n', s, '\n')



#a,s = align('GCTTCCTTTTCGCGAGCGT', 'AAGATTTTAGATT', 2.1393, 4.6378, -3.4524)
#print(a, '\n', s, '\n')

