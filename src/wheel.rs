use std::cmp;

use crate::constants::SEGMENT_LENGTH;
use crate::origin::Origin;
use crate::segment::WheelSegment;

/// Iterate through wheel_segment within [segment_start, segment_end)
///
/// When primes in this WheelSegment run out:
/// 1) If this reaches the overall end, then exit
/// 2) step in place to a new [segment_start, segment_end)
/// 3) create a new WheelSegment with these bounds and sieve it for primes.
pub struct Wheel {
    end: usize,
    origin_primes: Vec<usize>,
    segment_start: usize,
    segment_end: usize,
    wheel_segment: WheelSegment,
}

impl From<Origin> for Wheel {
    fn from(state: Origin) -> Wheel {
        let end = state.end;
        let origin_primes = state.origin_primes;

        let segment_start = cmp::min(state.start, end);
        let segment_end = cmp::min(segment_start + SEGMENT_LENGTH, end);
        let wheel_segment = WheelSegment::new(&origin_primes, segment_start, segment_end);

        Wheel {
            end,
            origin_primes,
            segment_start,
            segment_end,
            wheel_segment,
        }
    }
}

impl Iterator for Wheel {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.wheel_segment.next()
    }
}

impl Wheel {
    // Was this the last segment to sieve?
    pub fn done(&self) -> bool {
        self.segment_end == self.end
    }

    // Sieve the next segment, consuming and returning self
    pub fn advance_segment(mut self) -> Self {
        self.segment_start = self.segment_end;
        self.segment_end = cmp::min(self.segment_start + SEGMENT_LENGTH, self.end);
        self.wheel_segment =
            WheelSegment::new(&self.origin_primes, self.segment_start, self.segment_end);

        self
    }
}
