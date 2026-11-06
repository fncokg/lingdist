# Independently derived gold standards, reviewed and approved by the author.
# Based on the lingdist 3.0.0 calculation conventions reviewed on 2026-11-06.
# Never regenerate expected values from lingdist output.
# If the scientific definition changes, revise the derivation and constants
# together. These examples check public computational results only.
# Numeric comparisons allow an absolute error of 1e-10 per value; output row
# pairs and cost symbols are matched by name rather than incidental order.

# Default edit distance ----------------------------------------------------
# Split tokens at "_"; "pʰ" is one token. Matching costs 0; substitution,
# insertion and deletion each cost 1. Divide each word-pair distance by the
# longer token count, then average over jointly observed words, ignoring NA.
#
# Pair       word1    word2    word3     Mean
# L1 -> L2   1/3      1/3      missing   (1/3 + 1/3)/2 = 1/3
# L1 -> L3   1/3      1/2      missing   (1/3 + 1/2)/2 = 5/12
# L2 -> L3   2/3      1/3      1         (2/3 + 1/3 + 1)/3 = 2/3
# For L2 -> L3, word1 substitutes p with pʰ and deletes t; word2 deletes b.

golden_edit <- list(
    input = data.frame(
        word1 = c("pʰ_a_t", "p_a_t", "pʰ_a"),
        word2 = c("a_b", "a_b_c", "a_c"),
        word3 = c(NA_character_, "x", "y"),
        row.names = c("L1", "L2", "L3")
    ),
    expected = data.frame(
        lab1 = c("L1", "L1", "L2"), lab2 = c("L2", "L3", "L3"),
        dist = c(1/3, 5/12, 2/3)
    )
)

# Custom directed edit distance -------------------------------------------
# Matrix rows are source symbols and columns are targets; _NULL_ is the gap.
# a -> b: substitution costs 0.25, less than delete + insert = 0.75 + 0.6.
# b -> a: delete + insert costs 0.5 + 0.4 = 0.9, less than substitution = 1.5.
# Both self-distances are 0. All words have length 1, so normalization changes
# nothing. The b -> a value also checks the public string_edit_dist API.
golden_custom_edit <- list(
    input = data.frame(word = c("a", "b"), row.names = c("L1", "L2")),
    cost = as.data.frame(matrix(
        c(0, 0.25, 0.75, 1.5, 0, 0.5, 0.4, 0.6, 0),
        nrow = 3, byrow = TRUE,
        dimnames = list(c("a", "b", "_NULL_"), c("a", "b", "_NULL_"))
    )),
    expected = data.frame(
        lab1 = c("L1", "L1", "L2", "L2"),
        lab2 = c("L1", "L2", "L1", "L2"),
        dist = c(0, 1/4, 9/10, 0)
    )
)

# Weighted Jaccard distance ------------------------------------------------
# A single form's distance is the weight of the first differing category
# level. If one form is a prefix, use the first missing level's weight;
# identical forms have distance 0. Category weights are (1, 1/2, 1/4).
# Split forms at "#" and categories at "_". With form_strategy="weighting",
# pair forms by position, repeating the shorter side's last form. Normalize
# the used form weights (1, 1/2) to (2/3, 1/3).
#
# Pair       word1                           word2   word3     Mean
# L1 -> L2   (2/3)(1/4) + (1/3)(1) = 1/2      1/4     missing   3/8
# L1 -> L3   (2/3)(1)   + (1/3)(1) = 1        1/2     missing   3/4
# L2 -> L3   (2/3)(1)   + (1/3)(1/4) = 3/4    1/2     0         5/12
# C_1 vs C_1_z differs at the third level, giving 1/4. L1 and L3 have the
# same two forms in word1 but in opposite order; positional weighting gives 1.
golden_wjd <- list(
    input = data.frame(
        word1 = c("A_1_x#B_2", "A_1_y", "B_2#A_1_x"),
        word2 = c("C_1", "C_1_z", "C_2"),
        word3 = c(NA_character_, "D", "D"),
        row.names = c("L1", "L2", "L3")
    ),
    expected = data.frame(
        lab1 = c("L1", "L1", "L2"), lab2 = c("L2", "L3", "L3"),
        dist = c(3/8, 3/4, 5/12)
    )
)

# PMI: one update and convergence -----------------------------------------
# Both cases use the same two-row, four-word input. word2 and word3 are
# deliberately repeated to produce asymmetric counts and an intermediate
# learned cost. Split tokens at "_" with form_strategy="off".
# Set damping_alpha=1, smoothing_lambda=1, min_edit_cost=0.2, tol=1e-12 and
# normalize_method="longest". The smoothing and lower bound are chosen for
# readable arithmetic, not to reproduce the software's defaults.
# Initial matching costs are 0 and all nonmatching costs are 1. Below, g means
# _NULL_; counts refer to unordered nonmatching symbol pairs, including gaps.
#
# FIRST UPDATE
# ab -> ba has exactly three optimal alignments, each with total cost 2:
#
# Source   a b       a b g       g a b
# Target   b a       g b a       b a g
# Counts   ab:2      ag:2        bg:2
#
# Average equally over complete paths: each unordered pair receives 2/3.
# Matching edges affect path probabilities but add no PMI observations.
# The two a -> b words add 2 to ab. a -> aa has two optimal insertion positions;
# each contributes one ag, so its expected ag count is 1. Do not length-normalize
# training counts.
#
# Pair       Observed   After adding lambda=1
# ab         8/3        11/3
# ag         5/3         8/3
# bg         2/3         5/3
#
# Smoothed marginals are m(a)=19/3, m(b)=16/3, m(g)=13/3.
# Compute s(x,y)=ln(m(x)*m(y)/joint(x,y)):
#   s(ab)=ln(304/33) [minimum], s(ag)=ln(247/24), s(bg)=ln(208/15) [maximum].
# Probability denominators add a shared constant that cancels in the min-max
# mapping. Map these scores linearly to [0.2, 1], keeping the diagonal at 0:
#   cost(ab) = 0.2
#   cost(ag) = x = 0.2 + 0.8 * ln((247/24)/(304/33)) / ln((208/15)/(304/33))
#            = 0.41676891522749354
#   cost(bg) = 1
#
# Prediction uses the UPDATED costs. ab -> ba now prefers two substitutions
# (0.4) over two ag edits (2*x). The four length-normalized word distances are
# (0.2, 0.2, 0.2, x/2), with mean (0.6 + x/2)/4 = 0.20209611440343669.
# With max_epochs=1, converged=FALSE and a nonconvergence warning is expected.
#
# CONVERGED RESULT
# In the second training round, ab -> ba uniquely uses two substitutions.
# Observed counts are now ab=4, ag=1, bg=0. Add lambda to obtain (5, 2, 1),
# with marginals m(a)=7, m(b)=6, m(g)=3:
#   s(ab)=ln(42/5) [minimum], s(ag)=ln(21/2), s(bg)=ln(18) [maximum].
# After mapping to [0.2, 1]:
#   cost(ab) = 0.2
#   cost(ag) = y = 0.2 + 0.8 * ln(5/4)/ln(15/7) = 0.43422839486250131
#   cost(bg) = 1
# These costs preserve the optimal paths, so the third update changes nothing
# and meets the convergence tolerance. Allow max_epochs=20 and expect
# converged=TRUE without fixing the exact iteration count in the test.
# The four word distances are (0.2, 0.2, 0.2, y/2), with mean
# (0.6 + y/2)/4 = 0.20427854935781264.
#
# Both cases check the returned distance AND public 3x3 cost matrix: diagonal
# entries are 0, and off-diagonal entries are the symmetric ab, ag and bg costs.
golden_pmi <- list(
    input = data.frame(
        word1 = c("a_b", "b_a"), word2 = c("a", "b"),
        word3 = c("a", "b"), word4 = c("a", "a_a"),
        row.names = c("L1", "L2")
    ),
    # Explicit parameters keep this small example easy to audit.
    args = list(
        delim = "_", normalize_method = "longest", form_strategy = "off",
        damping_alpha = 1, smoothing_lambda = 1, min_edit_cost = 0.2,
        tol = 1e-12, detailed = FALSE, squareform = FALSE, quiet = TRUE
    ),
    one_epoch = list(
        expected = data.frame(lab1 = "L1", lab2 = "L2", dist = 0.20209611440343669),
        cost = matrix(
            c(0, 0.2, 0.41676891522749354,
              0.2, 0, 1,
              0.41676891522749354, 1, 0),
            nrow = 3, byrow = TRUE,
            dimnames = list(c("a", "b", "_NULL_"), c("a", "b", "_NULL_"))
        )
    ),
    converged = list(
        expected = data.frame(lab1 = "L1", lab2 = "L2", dist = 0.20427854935781264),
        cost = matrix(
            c(0, 0.2, 0.43422839486250131,
              0.2, 0, 1,
              0.43422839486250131, 1, 0),
            nrow = 3, byrow = TRUE,
            dimnames = list(c("a", "b", "_NULL_"), c("a", "b", "_NULL_"))
        )
    )
)
