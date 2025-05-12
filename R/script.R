#' Compute reverse complement of a DNA sequence
#'
#' @param sequence A single DNA sequence (character string), e.g., "ACGT"
#' @param complement_map Dictionary defining base complements (default: A:T, C:G)
#'
#' @return An object of class "ReverseComplementDNA" with attributes:
#' \itemize{
#'   \item \code{original_sequence}: the original input sequence
#'   \item \code{sequence_length}: sequence length
#'   \item \code{class}: "ReverseComplementDNA"
#' }
#' @examples
#' reverse_complement_dna("ACGT")
#' reverse_complement_dna("ACGTN", complement_map = c(A="T", T="A", C="G", G="C", N="N"))
reverse_complement_dna <- function(sequence, complement_map = c(A = "T", T = "A", C = "G", G = "C")) {
  stopifnot(is.character(sequence), length(sequence) == 1)

  chars <- unlist(strsplit(sequence, split = "")) # as vector because of mapping in complement_map[chars]

  if (!all(chars %in% names(complement_map))) {
    stop("Invalid characters in sequence, please check the complement_map.")
  }

  comp <- rev(complement_map[chars])
  result <- paste(comp, collapse = "") # vector to string

  class(result) <- "ReverseComplementDNA"
  attr(result, "original_sequence") <- sequence
  attr(result, "sequence_length") <- nchar(sequence)
  return(result)
}

#' One-hot encode a ReverseComplementDNA sequence
#'
#' @param rc_object An object of class "ReverseComplementDNA"
#' @param tokens A character vector of allowed tokens to encode (default: c("A", "C", "G", "T"))
#'
#' @return A matrix of size (length(tokens) x sequence_length) with one-hot encoding
#'
#' @examples
#' rc <- reverse_complement_dna("ACGT")
#' one_hot_encode_dna(rc)
one_hot_encode_dna <- function(rc_object, tokens = c("A", "C", "G", "T")) {
  stopifnot(inherits(rc_object, "ReverseComplementDNA"))
  stopifnot(is.character(tokens), length(tokens) > 0)

  seq <- as.character(rc_object)
  chars <- unlist(strsplit(seq, split = ""))
  if (!all(chars %in% tokens)) {
    stop("Sequence contains tokens not in the allowed token set.")
  }

  one_hot <- sapply(tokens, function(x) as.integer(chars == x))

  attr(one_hot, "original_class") <- class(rc_object)
  attr(one_hot, "original_sequence") <- attr(rc_object, "original_sequence")
  return(one_hot)
}

library("testthat")

# Tests for reverse_complement_dna()
test_that("reverse complement works for ACGT", {
  rc <- reverse_complement_dna("ACGT")
  expect_equal(as.character(rc), "ACGT")

})

test_that("reverse complement works for TTTG", {
  rc <- reverse_complement_dna("TTTG")
  expect_equal(as.character(rc), "CAAA")

})

test_that("reverse_complement_dna works correctly", {
  expect_equal(as.character(reverse_complement_dna("ACGT")), "ACGT")
  expect_equal(as.character(reverse_complement_dna("ACGTN", complement_map = c(A="T", T="A", C="G", G="C", N="N"))), "NACGT")
  expect_error(reverse_complement_dna("AXYZ")) # unknown characters
  expect_error(reverse_complement_dna(c("A", "T"))) # More than one sequence as an input
  expect_s3_class(reverse_complement_dna("ACGT"), "ReverseComplementDNA")
  expect_equal(attr(reverse_complement_dna("ACGT"), "original_sequence"), "ACGT")
})
rc <- reverse_complement_dna("ACGT")
encoded <- one_hot_encode_dna(rc)

test_that("one_hot_encode_dna works correctly", {
  rc <- reverse_complement_dna("ACGT")
  encoded <- one_hot_encode_dna(rc)

  expect_equal(dim(encoded), c(4, 4))

  expect_equal(unname(encoded[1, ]), c(1, 0, 0, 0))  # A
  expect_equal(unname(encoded[2, ]), c(0, 1, 0, 0))  # C
  expect_equal(unname(encoded[3, ]), c(0, 0, 1, 0))  # G
  expect_equal(unname(encoded[4, ]), c(0, 0, 0, 1))  # T

  expect_equal(attr(encoded, "original_class"), "ReverseComplementDNA")
  expect_equal(attr(encoded, "original_sequence"), "ACGT")
})

test_that("one_hot_encode_dna handles invalid input", {
  expect_error(one_hot_encode_dna("ACGT")) # is not S3 class
  rc <- reverse_complement_dna("ACGT")
  expect_error(one_hot_encode_dna(rc, tokens = c("X", "Y"))) # unknown tokens
})

#' Example DNA sequences
#'
#' A named character vector containing 10 randomly generated DNA sequences (length 20).
#'
#' @format Named character vector with 10 elements.
#' @usage data(example_sequences)
"example_sequences"

set.seed(123)
random_sequence <- function(n) {
  paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}

example_sequences <- setNames(
  replicate(10, random_sequence(20)),
  paste0("seq_", 1:10)
)

save(example_sequences, file = "example_sequences.RData")

#' Plot GC content for a set of DNA sequences
#'
#' @param sequences A named character vector of DNA sequences (names are IDs)
#'
#' @return A ggplot object showing GC content for each sequence
#' @import ggplot2
#' @examples
#' data(example_sequences)
#' plot_gc_content(example_sequences)
plot_gc_content <- function(sequences) {
  stopifnot(is.character(sequences), !is.null(names(sequences)))

  library(ggplot2)

  gc_content <- function(seq) {
    chars <- unlist(strsplit(seq, split = ""))
    sum(chars %in% c("G", "C")) / length(chars)
  }

  df <- data.frame(
    ID = names(sequences),
    GC_Content = sapply(sequences, gc_content)
  )

  p <- ggplot(df, aes(x = ID, y = GC_Content)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    ylim(0, 1) +
    theme_minimal() +
    labs(title = "GC Content per Sequence", x = "Sequence ID", y = "GC Content")

  return(p)
}

load("example_sequences.RData")
plot_gc_content(example_sequences)
rc = reverse_complement_dna(example_sequences[1])
one_hot_encode_dna(rc)
