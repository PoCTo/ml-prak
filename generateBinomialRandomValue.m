function x = generateBinomialRandomValue(n, p)
    x = sum(generateBernulliRandomSequence(p, n));