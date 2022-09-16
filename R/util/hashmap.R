library(digest)

# simple hashmap implementation that extends a demo package to account for 
# hash collisions.
# 
# source: https://github.com/nfultz/ht

#' Tiny Hash Table
#' 
#' This is a very basic implementation of a hash table using the \code{digest} package,
#' primarily for teaching functions and S3 for R programmers.
#' 
#' \code{ht} is an S3 class that extends \code{environment}, and additionally provides \code{[} and \code{[<-}.
#' It can use arbitrary R objects as keys and values.
#' 
#' 
#' @name ht
#' @seealso \code{\link[digest]{digest}}
#' @import digest
#' @examples
#'  x <- ht()
#'  x[1] <- 1
#'  x[1:2] <- 3:4
#'  x[1]
#'  x[1:2]
#'  1:2 %in.ht% x
#'  mget(ls(x),x)
#'  if(require(digest)) x[[digest(1:2)]]
NULL

#' Create ht structure
#' 
#' @export
#' @rdname ht
ht <- function(){
  structure(new.env(parent = emptyenv(), hash = TRUE), class="ht");
}

#' Export ht contents as a list of keys and values
#'
#' @param x an \code{ht} object
#' @export
#' @rdname ht
as.list.ht <- function(x) {
  # enumerate buckets in hashmap
  buckets = ls(x)
  # read all buckets
  contents = lapply(buckets, function(b) x[[b]])
  # merge key/value contents from across buckets
  list(
    keys = do.call(c, lapply(contents, function(x) x$keys)),
    values = do.call(c, lapply(contents, function(x) x$values))
  )
}

#' Retrieve value associated with key
#' 
#' @param x an \code{ht} object
#' @param key  A key object
#' @export
#' @rdname ht
`[.ht` <- function(x, key) {

  # initialize result with default value for missing keys
  res = NULL
  
  # bucket in which to search for key/value
  bucket = digest(key)
  
  if(exists(bucket, envir = x, inherits = FALSE)) {
    # retrieve keys/values in bucket
    contents = x[[bucket]]
    # search for key in bucket
    ind = which(sapply(contents$keys, function(k) identical(k, key)))
    # extract value for key
    if(length(ind) > 0) {
      res = contents$values[[ind]]
    }
  }
  
  res
}

#' Associate a value with a key
#' 
#' @param value A value object
#' @export
#' @rdname ht
`[<-.ht` <- function(x, key, value) {
  
  # bucket in which to store key/value
  bucket = digest(key)
  
  # build or update bucket for key/value
  if(exists(bucket, envir = x, inherits = FALSE)) {
    # retrieve keys/values in bucket
    contents = x[[bucket]]
    # search for key in bucket
    ind = which(sapply(contents$keys, function(k) identical(k, key)))
    # store key/value
    if(length(ind) > 0) {
      # update key's value in bucket
      contents$values[[ind]] = value
    } else {
      # append key/value pair to bucket (solution for hash collisions)
      contents$keys = c(contents$keys, list(key))
      contents$values = c(contents$values, list(value))
    }
  } else {
    # store key/value in new bucket
    contents = list(keys=list(key), values=list(value))
  }
  
  # update bucket contents in hashmap
  x[[bucket]] = contents
  
  x
}

#' Check for a key in a hashmap
#' 
#' @export
#' @rdname ht
`%in.ht%` <- function(key, x) {
  
  # initialize result with default value for missing keys
  res = FALSE
  
  # bucket in which to search for key/value
  bucket = digest(key)
  
  if(exists(bucket, envir = x, inherits = FALSE)) {
    # retrieve keys/values in bucket
    contents = x[[bucket]]
    # search for key in bucket
    ind = which(sapply(contents$keys, function(k) identical(k, key)))
    # extract value for key
    if(length(ind) > 0) {
      res = TRUE
    }
  }
  
  res
}