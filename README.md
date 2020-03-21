# Prime Generator

This crate generates prime numbers using the Sieve of Eratosthenes, with
several optimizations to take advantage of modern processors and memory. It's
by no means the fastest sieve (see Kim Walisch's
[primesieve](https://github.com/kimwalisch/primesieve/)), nor even
the fastest in Rust (see Huon Wilson's
[primal](https://docs.rs/primal/0.2.3/primal/)).

But it does fit my main purpose (solving Project Euler problems), and I think
this implementation is particularly legible and clean.

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.