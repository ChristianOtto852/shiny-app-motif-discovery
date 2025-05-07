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
            profile[j][i] = profile[j][i] / (t + 4)
    
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

def profile_randomly_generated_kmer(sequence, k, profile):
    probabilities = []
    for i in range(len(sequence) - k + 1):
        probability = 1
        for j in range(k):
            probability *= profile[j][sequence[i:(i+k)][j]]
        probabilities.append(probability)
    constant = sum(probabilities)
    new_probabilities = [x / constant for x in probabilities]
    start_motif = random.choices(range(len(sequence) - k + 1), weights = new_probabilities)[0]
    new_kmer = sequence[start_motif:(start_motif + k)]
    return new_kmer

def gibbs_sampler(dna: list[str], k: int, n: int) -> list[str]:
    """Implements the GibbsSampling algorithm for motif finding."""
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
            
        for j in range(n):
            i = random.randint(0, t-1)
            new_motif = motifs.pop(i)
            profile = profile_function(motifs, k)
            motif_i = profile_randomly_generated_kmer(dna[i], k, profile)
            motifs.insert(i, motif_i)
            if score_function(motifs, k) < score_function(best_motifs, k):
                best_motifs = motifs
    return best_motifs
