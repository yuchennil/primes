use std::cmp;
use std::mem;

/// Merge k min-sorted iterators I over values T into a single min-sorted iterator over T
///
/// This is equivalent to pushing the least element of each iterator into a priority queue and then
/// popping elements. But instead of using a BinaryHeap as the backing data store, we implement a
/// loser tournament tree. This lets us avoid the double O(log(k)) cost of a pop() followed by a
/// separate push().
///
/// The loser tournament tree has the following structure:
///                                     winner
///                                    loser[0]
///                loser[1]                                loser[2]
///      loser[3]            loser[4]            loser[5]            loser[6]
/// iters[0]  iters[1]  iters[2]  iters[3]  iters[4]
/// level_offset = 7
///
/// In this example tree, it's okay that iterators does not fill the last row. loser[6] exists and
/// will always be None, but that doesn't affect performance.
///
/// To iterate through the KWayMerge:
/// 1) Take winner's value
/// 2) pop the next value off winner's iterator
/// 3) starting from winner's iterator, replay tournament games in the brackets
/// 4) assign the overall tournament winner to the (now-vacant) winner
///
/// Since only winner's brackets have changed during this iteration, no other games need to be
/// replayed to maintain the heap property.
pub struct KWayMerge<T, I>
where
    T: Ord + Copy + Clone,
    I: Iterator<Item = T>,
{
    iterators: Vec<I>,
    level_offset: usize,
    winner: Node<T>,
    losers: Vec<Node<T>>,
}

impl<T, I> Iterator for KWayMerge<T, I>
where
    T: Ord + Copy + Clone,
    I: Iterator<Item = T>,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        let (value, index) = self.winner.0.take()?;
        let next_value = self.iterators[index].next();
        self.push(next_value, index);
        Some(value)
    }
}

impl<T, I> KWayMerge<T, I>
where
    T: Ord + Copy + Clone,
    I: Iterator<Item = T>,
{
    pub fn new(mut iterators: Vec<I>) -> KWayMerge<T, I> {
        let level_offset = KWayMerge::<T, I>::level_offset(iterators.len());
        let (winner, losers) = KWayMerge::initialize_tree(level_offset, &mut iterators);

        KWayMerge {
            iterators,
            level_offset,
            winner,
            losers,
        }
    }

    fn push(&mut self, option_value: Option<T>, index: usize) {
        let mut current_index = index + self.level_offset;
        let mut current_value = match option_value {
            Some(value) => Node(Some((value, index))),
            None => Node(None),
        };
        while let Some(parent_index) = KWayMerge::<T, I>::parent(current_index) {
            let parent_value = &mut self.losers[parent_index];
            if parent_value < &mut current_value {
                mem::swap(parent_value, &mut current_value);
            }
            current_index = parent_index;
        }
        // If winner had any existing value, it would be overwritten. self.push() should only
        // be called after taking self.winner.
        self.winner = current_value;
    }

    fn initialize_tree(level_offset: usize, iterators: &mut Vec<I>) -> (Node<T>, Vec<Node<T>>) {
        let winners_size = KWayMerge::<T, I>::left_child(level_offset);
        let mut winners = vec![Node(None); winners_size];
        for (index, spoke) in iterators.iter_mut().enumerate() {
            winners[index + level_offset] = Node(match spoke.next() {
                Some(value) => Some((value, index)),
                None => None,
            });
        }

        let mut losers = vec![Node(None); level_offset];
        let mut curr_offset = level_offset;
        while let Some(next_offset) = KWayMerge::<T, I>::parent(curr_offset) {
            for next_index in next_offset..curr_offset {
                let left_node = winners[KWayMerge::<T, I>::left_child(next_index)];
                let right_node = winners[KWayMerge::<T, I>::right_child(next_index)];
                let (winner, loser) = if left_node < right_node {
                    (left_node, right_node)
                } else {
                    (right_node, left_node)
                };
                if winner.0.is_some() {
                    winners[next_index] = winner;
                }
                if loser.0.is_some() {
                    losers[next_index] = loser;
                }
            }
            curr_offset = next_offset;
        }

        let winner = *winners.get(0).unwrap_or(&Node(None));

        (winner, losers)
    }

    fn level_offset(k: usize) -> usize {
        let mut power_of_two = 1;
        while power_of_two < k {
            power_of_two *= 2;
        }
        power_of_two - 1
    }
    fn parent(index: usize) -> Option<usize> {
        if index == 0 {
            return None;
        }
        Some((index - 1) / 2)
    }
    fn left_child(index: usize) -> usize {
        2 * index + 1
    }
    fn right_child(index: usize) -> usize {
        2 * index + 2
    }
}

#[derive(PartialEq, Copy, Clone)]
struct Node<T: Ord + Copy + Clone>(Option<(T, usize)>);

impl<T: Ord + Copy + Clone> PartialOrd for Node<T> {
    fn partial_cmp(&self, other: &Node<T>) -> Option<cmp::Ordering> {
        Some(match (self.0, other.0) {
            (Some(self_tuple), Some(other_tuple)) => self_tuple.cmp(&other_tuple),
            (Some(_), None) => cmp::Ordering::Less,
            (None, Some(_)) => cmp::Ordering::Greater,
            (None, None) => cmp::Ordering::Equal,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn k_way_merge_correct() {
        let sorted_array0 = vec![2, 3, 5, 7, 11];
        let sorted_array1 = vec![1, 1, 2, 3, 5, 8, 13];
        let iterators = vec![sorted_array0.into_iter(), sorted_array1.into_iter()];

        let k_way_merge = KWayMerge::new(iterators);
        assert_eq!(
            vec![1, 1, 2, 2, 3, 3, 5, 5, 7, 8, 11, 13],
            k_way_merge.collect::<Vec<_>>()
        );
    }
}
