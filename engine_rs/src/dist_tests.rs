    use super::*;

    // ── UniformBase ────────────────────────────────────────────────

    #[test]
    fn uniform_base_only_emits_canonical_bases() {
        let mut rng = Rng::new(1);
        let dist = UniformBase;
        for _ in 0..1000 {
            let b = dist.sample(&mut rng);
            assert!(
                matches!(b, b'A' | b'C' | b'G' | b'T'),
                "UniformBase emitted unexpected byte 0x{:02x}",
                b
            );
        }
    }

    #[test]
    fn uniform_base_covers_all_four_bases() {
        let mut rng = Rng::new(123);
        let dist = UniformBase;
        let mut seen = [false; 4];
        for _ in 0..1000 {
            let b = dist.sample(&mut rng);
            let idx = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("unexpected base"),
            };
            seen[idx] = true;
        }
        for (i, &was_seen) in seen.iter().enumerate() {
            assert!(was_seen, "base index {} never appeared in 1000 draws", i);
        }
    }

    #[test]
    fn uniform_base_distribution_is_roughly_uniform() {
        let mut rng = Rng::new(0xfeed);
        let dist = UniformBase;
        let mut counts = [0u32; 4];
        let n = 10_000u32;
        for _ in 0..n {
            let idx = match dist.sample(&mut rng) {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => unreachable!(),
            };
            counts[idx] += 1;
        }
        // Expected count per bucket = 2500. Generous tolerance.
        for &c in &counts {
            assert!(
                (2200..=2800).contains(&c),
                "UniformBase bucket count {} outside [2200, 2800]",
                c
            );
        }
    }

    #[test]
    fn uniform_base_same_seed_same_stream() {
        let mut a = Rng::new(7);
        let mut b = Rng::new(7);
        let dist = UniformBase;
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
        }
    }

    // ── UniformInt ─────────────────────────────────────────────────

    #[test]
    #[should_panic(expected = "UniformInt: max")]
    fn uniform_int_new_rejects_max_le_min() {
        let _ = UniformInt::new(5, 5);
    }

    #[test]
    #[should_panic(expected = "UniformInt: max")]
    fn uniform_int_new_rejects_inverted_range() {
        let _ = UniformInt::new(10, 5);
    }

    #[test]
    #[should_panic(expected = "UniformInt: span")]
    fn uniform_int_new_rejects_oversized_span() {
        // Span exceeds u32::MAX → must panic explicitly, not silently
        // truncate. Pre-cleanup the constructor accepted this and
        // sample() biased the result; now it's caught at construction.
        let _ = UniformInt::new(0, (u32::MAX as i64) + 2);
    }

    #[test]
    fn uniform_int_new_accepts_max_span() {
        // Right at the boundary — span == u32::MAX is allowed.
        let dist = UniformInt::new(0, u32::MAX as i64);
        assert_eq!(dist.span(), u32::MAX as u64);
    }

    #[test]
    fn uniform_int_stays_in_bounds() {
        let mut rng = Rng::new(42);
        let dist = UniformInt::new(3, 10); // [3, 10) → 3..=9
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            assert!(v >= 3, "UniformInt(3,10) produced {} (< 3)", v);
            assert!(v < 10, "UniformInt(3,10) produced {} (>= 10)", v);
        }
    }

    #[test]
    fn uniform_int_covers_full_range() {
        let mut rng = Rng::new(99);
        let dist = UniformInt::new(0, 5);
        let mut seen = [false; 5];
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            seen[v as usize] = true;
        }
        for (i, &was_seen) in seen.iter().enumerate() {
            assert!(was_seen, "UniformInt(0,5) never produced {}", i);
        }
    }

    #[test]
    fn uniform_int_negative_min_works() {
        let mut rng = Rng::new(11);
        let dist = UniformInt::new(-3, 4); // [-3, 4) → -3, -2, -1, 0, 1, 2, 3
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            assert!(v >= -3, "{}", v);
            assert!(v < 4, "{}", v);
        }
    }

    #[test]
    fn uniform_int_span_one_returns_min() {
        let mut rng = Rng::new(1);
        let dist = UniformInt::new(7, 8);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut rng), 7);
        }
    }

    #[test]
    fn uniform_int_same_seed_same_stream() {
        let mut a = Rng::new(0xbabe);
        let mut b = Rng::new(0xbabe);
        let dist = UniformInt::new(0, 100);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
        }
    }

    #[test]
    fn uniform_int_accessors_round_trip() {
        let dist = UniformInt::new(-5, 12);
        assert_eq!(dist.min(), -5);
        assert_eq!(dist.max(), 12);
        assert_eq!(dist.span(), 17);
    }

    // ── Trait object usage ─────────────────────────────────────────

    #[test]
    fn box_dyn_distribution_dispatches_correctly() {
        // The whole point of D4: heterogeneous distributions stored
        // as trait objects. This compiles only if the trait is
        // actually dyn-compatible with the chosen Output type.
        let mut rng = Rng::new(1);

        let base_dist: Box<dyn Distribution<Output = u8>> = Box::new(UniformBase);
        let int_dist: Box<dyn Distribution<Output = i64>> = Box::new(UniformInt::new(0, 10));

        let b = base_dist.sample(&mut rng);
        let i = int_dist.sample(&mut rng);

        assert!(matches!(b, b'A' | b'C' | b'G' | b'T'));
        assert!((0..10).contains(&i));
    }

    #[test]
    fn vec_of_homogeneous_dyn_distributions() {
        // Storing many distributions in a heterogeneous container.
        let dists: Vec<Box<dyn Distribution<Output = u8>>> = vec![
            Box::new(UniformBase),
            Box::new(UniformBase),
            Box::new(UniformBase),
        ];
        let mut rng = Rng::new(0);
        for d in &dists {
            let b = d.sample(&mut rng);
            assert!(matches!(b, b'A' | b'C' | b'G' | b'T'));
        }
    }

    #[derive(Clone, Debug)]
    struct NoSupportDist;

    impl Distribution for NoSupportDist {
        type Output = i64;

        fn sample(&self, _rng: &mut Rng) -> i64 {
            0
        }
    }

    #[test]
    fn sample_filtered_result_reports_support_unavailable() {
        let mut rng = Rng::new(1);
        let err = sample_filtered_result(&mut rng, &NoSupportDist, |_| true).unwrap_err();
        assert_eq!(err, FilteredSampleError::SupportUnavailable);
    }

    #[test]
    fn sample_filtered_result_reports_empty_admissible_support() {
        let mut rng = Rng::new(1);
        let dist = UniformInt::new(0, 4);
        let err = sample_filtered_result(&mut rng, &dist, |_| false).unwrap_err();
        assert_eq!(err, FilteredSampleError::EmptyAdmissibleSupport);
    }

    #[test]
    fn sample_filtered_result_samples_from_admissible_subset() {
        let mut rng = Rng::new(1);
        let dist = UniformInt::new(0, 10);
        for _ in 0..100 {
            let value = sample_filtered_result(&mut rng, &dist, |v| *v >= 7).unwrap();
            assert!((7..10).contains(&value));
        }
    }

    #[test]
    fn sample_filtered_permissive_collapses_filter_errors_to_none() {
        let mut rng = Rng::new(1);
        let dist = UniformInt::new(0, 4);
        assert_eq!(sample_filtered(&mut rng, &dist, |_| false), None);
        assert_eq!(sample_filtered(&mut rng, &NoSupportDist, |_| true), None);
    }

    // ── EmpiricalLengthDist ────────────────────────────────────────

    #[test]
    #[should_panic(expected = "must have at least one entry")]
    fn empirical_rejects_empty_histogram() {
        let _ = EmpiricalLengthDist::from_pairs(std::iter::empty());
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn empirical_rejects_zero_weight() {
        let _ = EmpiricalLengthDist::from_pairs(vec![(5, 0.0), (6, 1.0)]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn empirical_rejects_negative_weight() {
        let _ = EmpiricalLengthDist::from_pairs(vec![(5, -1.0)]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn empirical_rejects_nan_weight() {
        let _ = EmpiricalLengthDist::from_pairs(vec![(5, f64::NAN)]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn empirical_rejects_infinite_weight() {
        let _ = EmpiricalLengthDist::from_pairs(vec![(5, f64::INFINITY)]);
    }

    #[test]
    fn empirical_single_value_always_returned() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(42, 1.0)]);
        let mut rng = Rng::new(0xc0ff_ee);
        for _ in 0..1000 {
            assert_eq!(dist.sample(&mut rng), 42);
        }
    }

    #[test]
    fn empirical_construction_accessors_round_trip() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 2.0), (2, 1.0)]);
        assert_eq!(dist.len(), 3);
        assert!(!dist.is_empty());
        assert_eq!(dist.values(), &[0, 1, 2]);
        assert!((dist.total_weight() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn empirical_from_values_and_weights_matches_from_pairs() {
        let pairs = EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 2.0), (2, 3.0)]);
        let split =
            EmpiricalLengthDist::from_values_and_weights(vec![0, 1, 2], vec![1.0, 2.0, 3.0]);

        let mut a = Rng::new(7);
        let mut b = Rng::new(7);
        for _ in 0..100 {
            assert_eq!(pairs.sample(&mut a), split.sample(&mut b));
        }
    }

    #[test]
    #[should_panic(expected = "values and weights must have equal lengths")]
    fn empirical_from_values_and_weights_rejects_length_mismatch() {
        let _ = EmpiricalLengthDist::from_values_and_weights(vec![0, 1], vec![1.0]);
    }

    #[test]
    fn empirical_stays_in_value_set() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(3, 1.0), (7, 2.0), (11, 1.0), (-2, 0.5)]);
        let mut rng = Rng::new(17);
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            assert!(matches!(v, 3 | 7 | 11 | -2));
        }
    }

    #[test]
    fn empirical_covers_full_value_set() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0)]);
        let mut rng = Rng::new(99);
        let mut seen = [false; 4];
        for _ in 0..1000 {
            let v = dist.sample(&mut rng);
            seen[v as usize] = true;
        }
        for (i, &was_seen) in seen.iter().enumerate() {
            assert!(was_seen, "value {} never appeared", i);
        }
    }

    #[test]
    fn empirical_weights_are_respected() {
        // 90/10 weight split — over 10,000 draws the heavy value
        // should appear roughly 9,000 times. Generous tolerance
        // catches statistical jitter.
        let dist = EmpiricalLengthDist::from_pairs(vec![(0, 0.9), (1, 0.1)]);
        let mut rng = Rng::new(0xabba);
        let mut zero_count = 0;
        let n = 10_000;
        for _ in 0..n {
            if dist.sample(&mut rng) == 0 {
                zero_count += 1;
            }
        }
        assert!(
            (8500..=9500).contains(&zero_count),
            "expected zero_count ~9000 of 10000, got {}",
            zero_count
        );
    }

    #[test]
    fn empirical_uniform_weights_produce_roughly_uniform_distribution() {
        let dist =
            EmpiricalLengthDist::from_pairs(vec![(0, 1.0), (1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0)]);
        let mut rng = Rng::new(0xdead);
        let mut counts = [0u32; 5];
        let n = 10_000u32;
        for _ in 0..n {
            counts[dist.sample(&mut rng) as usize] += 1;
        }
        // Expected ~2000 per bucket. Generous tolerance.
        for (i, &c) in counts.iter().enumerate() {
            assert!(
                (1700..=2300).contains(&c),
                "bucket {} count {} outside [1700, 2300]",
                i,
                c
            );
        }
    }

    #[test]
    fn empirical_same_seed_same_stream() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(10, 1.0), (20, 2.0), (30, 3.0)]);
        let mut a = Rng::new(0xfeed);
        let mut b = Rng::new(0xfeed);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
        }
    }

    #[test]
    fn empirical_negative_values_supported() {
        let dist = EmpiricalLengthDist::from_pairs(vec![(-5, 1.0), (0, 1.0), (5, 1.0)]);
        let mut rng = Rng::new(13);
        for _ in 0..100 {
            let v = dist.sample(&mut rng);
            assert!(matches!(v, -5 | 0 | 5));
        }
    }

    #[test]
    fn empirical_works_through_box_dyn() {
        let dist: Box<dyn Distribution<Output = i64>> =
            Box::new(EmpiricalLengthDist::from_pairs(vec![(7, 1.0)]));
        let mut rng = Rng::new(0);
        assert_eq!(dist.sample(&mut rng), 7);
    }

    // ── AllelePoolDist ─────────────────────────────────────────────

    use crate::ir::Segment;
    use crate::refdata::Allele;

    /// Build a pool of `n` named alleles for testing. The base
    /// sequences are all single-byte `b'A'` since we don't care
    /// about content here, only sampling behavior.
    fn make_pool(n: usize) -> AllelePool {
        let mut p = AllelePool::new();
        for i in 0..n {
            let _ = p.push(Allele {
                name: format!("a{}*01", i),
                gene: format!("a{}", i),
                seq: b"A".to_vec(),
                segment: Segment::V,
                anchor: None,
            });
        }
        p
    }

    #[test]
    #[should_panic(expected = "pool must contain at least one allele")]
    fn allele_pool_dist_uniform_rejects_empty_pool() {
        let p = AllelePool::new();
        let _ = AllelePoolDist::uniform(&p);
    }

    #[test]
    #[should_panic(expected = "weights")]
    fn allele_pool_dist_from_weights_rejects_size_mismatch() {
        let p = make_pool(3);
        let _ = AllelePoolDist::from_weights(&p, vec![1.0, 2.0]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn allele_pool_dist_from_weights_rejects_zero_weight() {
        let p = make_pool(2);
        let _ = AllelePoolDist::from_weights(&p, vec![0.0, 1.0]);
    }

    #[test]
    #[should_panic(expected = "weight at index")]
    fn allele_pool_dist_from_weights_rejects_nan() {
        let p = make_pool(1);
        let _ = AllelePoolDist::from_weights(&p, vec![f64::NAN]);
    }

    #[test]
    fn allele_pool_dist_uniform_construction_round_trip() {
        let p = make_pool(5);
        let dist = AllelePoolDist::uniform(&p);
        assert_eq!(dist.len(), 5);
        assert!(!dist.is_empty());
        assert!((dist.total_weight() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn allele_pool_dist_uniform_covers_all_alleles() {
        let p = make_pool(10);
        let dist = AllelePoolDist::uniform(&p);
        let mut rng = Rng::new(99);
        let mut seen = vec![false; 10];
        for _ in 0..2000 {
            let id = dist.sample(&mut rng);
            assert!(id.as_usize() < 10, "out-of-bounds AlleleId");
            seen[id.as_usize()] = true;
        }
        for (i, &was_seen) in seen.iter().enumerate() {
            assert!(was_seen, "AlleleId({}) never sampled", i);
        }
    }

    #[test]
    fn allele_pool_dist_uniform_is_roughly_uniform() {
        let p = make_pool(4);
        let dist = AllelePoolDist::uniform(&p);
        let mut rng = Rng::new(0xfade);
        let mut counts = [0u32; 4];
        let n = 10_000;
        for _ in 0..n {
            counts[dist.sample(&mut rng).as_usize()] += 1;
        }
        // Expected ~2500 per bucket.
        for (i, &c) in counts.iter().enumerate() {
            assert!(
                (2200..=2800).contains(&c),
                "bucket {} count {} outside [2200, 2800]",
                i,
                c
            );
        }
    }

    #[test]
    fn allele_pool_dist_weights_are_respected() {
        // 80/10/10 weighted pool. Heavy allele should dominate.
        let p = make_pool(3);
        let dist = AllelePoolDist::from_weights(&p, vec![0.8, 0.1, 0.1]);
        let mut rng = Rng::new(0xfeed);
        let mut counts = [0u32; 3];
        let n = 10_000;
        for _ in 0..n {
            counts[dist.sample(&mut rng).as_usize()] += 1;
        }
        // Allele 0 ~80% (8000); alleles 1 and 2 ~10% each (1000).
        assert!(
            (7500..=8500).contains(&counts[0]),
            "heavy bucket count {} outside [7500, 8500]",
            counts[0]
        );
        assert!(
            (700..=1300).contains(&counts[1]) && (700..=1300).contains(&counts[2]),
            "light bucket counts {:?} outside [700, 1300]",
            &counts[1..]
        );
    }

    #[test]
    fn allele_pool_dist_same_seed_same_stream() {
        let p = make_pool(7);
        let dist = AllelePoolDist::uniform(&p);
        let mut a = Rng::new(0xcafe);
        let mut b = Rng::new(0xcafe);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut a), dist.sample(&mut b));
        }
    }

    #[test]
    fn allele_pool_dist_single_allele_always_returned() {
        let p = make_pool(1);
        let dist = AllelePoolDist::uniform(&p);
        let mut rng = Rng::new(1);
        for _ in 0..100 {
            assert_eq!(dist.sample(&mut rng), AlleleId::new(0));
        }
    }

    #[test]
    fn allele_pool_dist_works_through_box_dyn() {
        let p = make_pool(3);
        let dist: Box<dyn Distribution<Output = AlleleId>> = Box::new(AllelePoolDist::uniform(&p));
        let mut rng = Rng::new(0);
        let id = dist.sample(&mut rng);
        assert!(id.as_usize() < 3);
    }

    #[test]
    fn allele_pool_dist_sampled_ids_resolve_in_pool() {
        // Integration sanity: sampled AlleleIds round-trip through
        // the pool. This is the structural guarantee the dist gives:
        // sample → get → Some(allele).
        let p = make_pool(5);
        let dist = AllelePoolDist::uniform(&p);
        let mut rng = Rng::new(42);
        for _ in 0..50 {
            let id = dist.sample(&mut rng);
            let resolved = p.get(id);
            assert!(resolved.is_some(), "sampled id {:?} did not resolve", id);
        }
    }

    // ── restricted_uniform tests ────────────────────────────────────

    #[test]
    fn allele_pool_dist_restricted_uniform_only_samples_allowed_ids() {
        let p = make_pool(10);
        let allowed = vec![AlleleId::new(2), AlleleId::new(5), AlleleId::new(7)];
        let dist = AllelePoolDist::restricted_uniform(&p, allowed.clone());
        assert_eq!(dist.len(), 3);
        let mut rng = Rng::new(0xbeef);
        for _ in 0..500 {
            let id = dist.sample(&mut rng);
            assert!(
                allowed.contains(&id),
                "sampled id {:?} not in allowed set",
                id
            );
        }
    }

    #[test]
    fn allele_pool_dist_restricted_uniform_covers_all_allowed() {
        let p = make_pool(8);
        let allowed = vec![AlleleId::new(1), AlleleId::new(4), AlleleId::new(6)];
        let dist = AllelePoolDist::restricted_uniform(&p, allowed.clone());
        let mut rng = Rng::new(0xc0de);
        let mut seen = vec![false; allowed.len()];
        for _ in 0..2000 {
            let id = dist.sample(&mut rng);
            let pos = allowed.iter().position(|&a| a == id).unwrap();
            seen[pos] = true;
        }
        assert!(seen.iter().all(|&s| s), "not every allowed id was sampled");
    }

    #[test]
    fn allele_pool_dist_restricted_uniform_single_id_always_returns_it() {
        let p = make_pool(20);
        let dist = AllelePoolDist::restricted_uniform(&p, vec![AlleleId::new(13)]);
        let mut rng = Rng::new(7);
        for _ in 0..50 {
            assert_eq!(dist.sample(&mut rng), AlleleId::new(13));
        }
    }

    #[test]
    fn allele_pool_dist_restricted_uniform_support_round_trip() {
        let p = make_pool(6);
        let allowed = vec![AlleleId::new(0), AlleleId::new(3), AlleleId::new(5)];
        let dist = AllelePoolDist::restricted_uniform(&p, allowed.clone());
        let support = dist.support().expect("support is Some");
        let support_ids: Vec<AlleleId> = support.iter().map(|(id, _)| *id).collect();
        assert_eq!(support_ids, allowed);
        for (_, weight) in &support {
            assert!((*weight - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    #[should_panic(expected = "allowed_ids must be non-empty")]
    fn allele_pool_dist_restricted_uniform_rejects_empty() {
        let p = make_pool(3);
        let _ = AllelePoolDist::restricted_uniform(&p, vec![]);
    }

    #[test]
    #[should_panic(expected = "out of range")]
    fn allele_pool_dist_restricted_uniform_rejects_out_of_range_id() {
        let p = make_pool(3);
        let _ = AllelePoolDist::restricted_uniform(&p, vec![AlleleId::new(5)]);
    }

    #[test]
    #[should_panic(expected = "duplicate")]
    fn allele_pool_dist_restricted_uniform_rejects_duplicate_ids() {
        let p = make_pool(3);
        let _ = AllelePoolDist::restricted_uniform(&p, vec![AlleleId::new(1), AlleleId::new(1)]);
    }
