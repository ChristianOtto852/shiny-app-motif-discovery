def profile_most_probable_kmer(text: str, k: int,
                               profile: list[dict[str, float]]) -> str:
    """Identifies the most probable k-mer according to a given profile matrix.

    The profile matrix is represented as a list of columns, where the i-th element is a map
    whose keys are strings ("A", "C", "G", and "T") and whose values represent the probability
    associated with this symbol in the i-th column of the profile matrix.
    """
    profile_most_probable_index = 0
    profile_probability = 0
    for i in range(len(text) - k + 1):
        probability = 1
        for j in range(k):
            probability *= profile[j][text[i:(i+k)][j]]
        if probability > profile_probability:
            profile_probability = probability
            profile_most_probable_index = i
    
    return text[profile_most_probable_index:(profile_most_probable_index + k)]

def profile_function(motifs, k):
    t = len(motifs)
    profile = [[1] * k for i in range(4)]
    
    for i in range(k):
        for j in range(t):
            if motifs[j][i] == "A":
                profile[0][i] += 1
            if motifs[j][i] == "C":
                profile[1][i] += 1
            if motifs[j][i] == "G":
                profile[2][i] += 1
            if motifs[j][i] == "T":
                profile[3][i] += 1
    for i in range(k):
        for j in range(4):
            profile[j][i] = profile[j][i] / (t+4)
    
    key = ["A", "C", "G", "T"]
    profile_dictionary = []
    for i in range(k):
        temp = {}
        for j in range(4):
            temp[key[j]] = profile[j][i]
        profile_dictionary.append(temp)
    return profile_dictionary

def score_function(motifs, k):
    t = len(motifs)
    profile = [[0] * k for i in range(4)]
    
    for i in range(k):
        for j in range(t):
            if motifs[j][i] == "A":
                profile[0][i] += 1
            if motifs[j][i] == "C":
                profile[1][i] += 1
            if motifs[j][i] == "G":
                profile[2][i] += 1
            if motifs[j][i] == "T":
                profile[3][i] += 1
    
    score = 0
    for i in range(k):
        Max = 0
        for j in range(4):
            if profile[j][i] > Max:
                Max = profile[j][i]
        score += t - Max
    return score

def greedy_motif_search_pseudocounts(dna: list[str], k: int) -> list[str]:
    """Augments the GreedyMotifSearch algorithm with pseudocounts."""
    t = len(dna)
    best_motifs = [strings[0:k] for strings in dna]
    
    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i:(i+k)]]
        for j in range(1,t):
            profile = profile_function(motifs, k)
            motif = profile_most_probable_kmer(dna[j], k, profile)
            motifs.append(motif)
        if score_function(motifs, k) < score_function(best_motifs, k):
            best_motifs = motifs
    return best_motifs
