source(test_path("fixtures.R"), local = TRUE)

expect_golden_numbers <- function(actual, expected) {
    expect_identical(dim(actual), dim(expected))
    expect_length(actual, length(expected))
    expect_true(all(is.finite(actual)))
    expect_lte(max(abs(actual - expected)), 1e-10)
}

expect_golden_distances <- function(actual, expected) {
    # Only the labeled row pairs matter, not their return order.
    key <- function(x) paste(x$lab1, x$lab2, sep = "/")
    actual_keys <- key(actual)
    expected_keys <- key(expected)
    expect_identical(anyDuplicated(actual_keys), 0L)
    expect_length(actual_keys, length(expected_keys))
    expect_setequal(actual_keys, expected_keys)
    expect_golden_numbers(actual$dist[match(expected_keys, actual_keys)], expected$dist)
}

expect_golden_pmi <- function(actual, expected) {
    expect_golden_distances(actual$result, expected$expected)
    symbols <- rownames(expected$cost)
    expect_setequal(rownames(actual$cost), symbols)
    expect_setequal(colnames(actual$cost), symbols)
    expect_golden_numbers(as.matrix(actual$cost[symbols, symbols]), expected$cost)
}

test_that("default edit distances match the reviewed gold standard", {
    for (n_threads in c(1L, 2L)) {
        actual <- pw_edit_dist(
            golden_edit$input, delim = "_", normalize_method = "longest",
            form_strategy = "off", symmetric = TRUE, detailed = FALSE,
            squareform = FALSE, default_sub_cost = 1, default_ins_del_cost = 1,
            n_threads = n_threads, quiet = TRUE
        )
        expect_golden_distances(actual, golden_edit$expected)
    }
})

test_that("custom edit costs match the reviewed directed distances", {
    for (n_threads in c(1L, 2L)) {
        actual <- pw_edit_dist(
            golden_custom_edit$input, cost_mat = golden_custom_edit$cost,
            delim = "_", normalize_method = "longest", form_strategy = "off",
            symmetric = FALSE, detailed = FALSE, squareform = FALSE,
            n_threads = n_threads, quiet = TRUE
        )
        expect_golden_distances(actual, golden_custom_edit$expected)
    }
    actual <- string_edit_dist(
        "b", "a", cost_mat = golden_custom_edit$cost, delim = "_",
        normalize_method = "longest", return_alignments = FALSE, quiet = TRUE
    )
    expect_golden_numbers(actual$distance, 9/10)
})

test_that("WJD matches the reviewed gold standard", {
    for (n_threads in c(1L, 2L)) {
        actual <- pw_wjd(
            golden_wjd$input, cate_delim = "_", cate_weights = c(1, 1/2, 1/4),
            form_strategy = "weighting", form_delim = "#", form_weights = c(1, 1/2),
            detailed = FALSE, squareform = FALSE, n_threads = n_threads, quiet = TRUE
        )
        expect_golden_distances(actual, golden_wjd$expected)
    }
})

test_that("one PMI update matches the reviewed distances and learned costs", {
    for (n_threads in c(1L, 2L)) {
        expect_warning(
            actual <- do.call(pw_pmi_dist, c(
                list(data = golden_pmi$input), golden_pmi$args,
                list(max_epochs = 1L, n_threads = n_threads)
            )),
            "PMI distance computation did not converge within the maximum number of epochs.",
            fixed = TRUE
        )
        expect_golden_pmi(actual, golden_pmi$one_epoch)
        expect_false(actual$converged)
    }
})

test_that("PMI converges to the reviewed distances and learned costs", {
    for (n_threads in c(1L, 2L)) {
        actual <- do.call(pw_pmi_dist, c(
            list(data = golden_pmi$input), golden_pmi$args,
            list(max_epochs = 20L, n_threads = n_threads)
        ))
        expect_golden_pmi(actual, golden_pmi$converged)
        expect_true(actual$converged)
    }
})
