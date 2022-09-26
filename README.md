# Range Minimum Query (RMQ) [![Crates.io][crates-badge]][crates-url] [![Docs.rs][docs-badge]][docs-rs] [![MIT licensed][mit-badge]][mit-url]

[crates-badge]: https://img.shields.io/crates/v/range_minimum_query.svg
[crates-url]: https://crates.io/crates/range_minimum_query
[mit-badge]: https://img.shields.io/badge/license-MIT-blue.svg
[mit-url]: https://opensource.org/licenses/MIT
[docs-rs]: https://docs.rs/range_minimum_query
[docs-badge]: https://img.shields.io/docsrs/range_minimum_query/0.1.0

This crate implements a succinct data structure to solve the Range Minimum Query (RMQ) problem in constant time and linear space.

This code was derived (almost verbatim) from a succinct version of RMQ, implemented by Giuseppe Ottaviano's succinct c++ library: https://github.com/ot/succinct.

# What is RMQ (taken from Wikipedia)

In computer science, a range minimum query (RMQ) solves the problem of finding the minimal value in a sub-array of an array of comparable objects. Range minimum queries have several use cases in computer science, such as the lowest common ancestor problem and the longest common prefix problem (LCP). 

For example, when `A = [0,5,2,5,4,3,1,6,3]`, then the answer to the range minimum query for the sub-array `A[2...7] = [2,5,4,3,1,6]` is `6`, as `A[6] = 1`. 

# Usage

```rust
use range_minimum_query::Rmq;

let a = [0,5,2,5,4,3,1,6,3];
let rmq = Rmq::from_iter(a);
let res = rmq.range_minimum(2..=7);
assert_eq!(res.unwrap(),6);
```

# License

MIT
