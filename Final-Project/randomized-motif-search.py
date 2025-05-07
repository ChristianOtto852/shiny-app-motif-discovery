import random

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

def randomized_motif_search(dna: list[str], k: int) -> list[str]:
    """Implements the RandomizedMotifSearch algorithm with pseudocounts."""
    motifs = []
    t = len(dna)
    for string in dna:
        kmer_start = random.randint(0, len(string) - k)
        motifs.append(string[kmer_start:(kmer_start + k)])
    best_motifs = motifs
    
    for sims in range(1000):
        motifs = []
        for string in dna:
            kmer_start = random.randint(0, len(string) - k)
            motifs.append(string[kmer_start:(kmer_start + k)])

        while True:
            profile = profile_function(motifs, k)
            motifs = []
            for string in dna:
                motifs.append(profile_most_probable_kmer(string, k, profile))
            if score_function(motifs, k) < score_function(best_motifs, k):
                best_motifs = motifs
            else:
                break
    
    return best_motifs
