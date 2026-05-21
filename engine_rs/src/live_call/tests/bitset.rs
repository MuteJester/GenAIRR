use super::super::AlleleBitSet;
use super::id;

#[test]
fn allele_bitset_insert_remove_and_iterate_are_stable() {
    let mut set = AlleleBitSet::empty(130);
    assert!(set.is_empty());

    assert!(set.insert(id(0)));
    assert!(set.insert(id(64)));
    assert!(set.insert(id(129)));
    assert!(!set.insert(id(64)));

    assert_eq!(set.len(), 3);
    assert!(set.contains(id(0)));
    assert!(set.contains(id(64)));
    assert!(set.contains(id(129)));
    assert_eq!(set.to_ids(), vec![id(0), id(64), id(129)]);

    assert!(set.remove(id(64)));
    assert!(!set.remove(id(64)));
    assert_eq!(set.to_ids(), vec![id(0), id(129)]);
}

#[test]
fn allele_bitset_union_and_intersection_require_same_universe() {
    let a = AlleleBitSet::from_ids(8, [id(1), id(2), id(5)]);
    let b = AlleleBitSet::from_ids(8, [id(2), id(3), id(5)]);

    assert_eq!(a.unioned(&b).to_ids(), vec![id(1), id(2), id(3), id(5)]);
    assert_eq!(a.intersected(&b).to_ids(), vec![id(2), id(5)]);
}

#[test]
#[should_panic(expected = "outside universe length")]
fn allele_bitset_rejects_out_of_universe_id() {
    let mut set = AlleleBitSet::empty(2);
    set.insert(id(2));
}

#[test]
fn full_bitset_masks_unused_tail_bits() {
    let set = AlleleBitSet::full(65);
    assert_eq!(set.len(), 65);
    assert_eq!(set.to_ids().first(), Some(&id(0)));
    assert_eq!(set.to_ids().last(), Some(&id(64)));
}
