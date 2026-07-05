    use super::*;

    fn pool_with(genes: &[(&str, &str)]) -> AllelePool {
        // genes: (allele_name, gene_name)
        let mut p = AllelePool::new();
        for (name, gene) in genes {
            let _ = p.push(Allele {
                name: (*name).to_string(),
                gene: (*gene).to_string(),
                seq: vec![b'A'; 10],
                segment: Segment::V,
                anchor: Some(3),
                functional_status: None,
                subregions: Vec::new(),
            });
        }
        p
    }

    #[test]
    fn gene_index_groups_alleles_by_gene_in_first_appearance_order() {
        let pool = pool_with(&[
            ("IGHV1-2*01", "IGHV1-2"),
            ("IGHV1-2*02", "IGHV1-2"),
            ("IGHV3-23*01", "IGHV3-23"),
        ]);
        let idx = GeneIndex::build(&pool);

        assert_eq!(idx.len(), 2);
        let g12 = idx.gene_id("IGHV1-2").unwrap();
        let g323 = idx.gene_id("IGHV3-23").unwrap();
        assert_eq!(g12.index(), 0); // first appearance
        assert_eq!(g323.index(), 1);
        assert_eq!(idx.gene_name(g12), "IGHV1-2");
        assert_eq!(idx.alleles_of(g12), &[AlleleId::new(0), AlleleId::new(1)]);
        assert_eq!(idx.alleles_of(g323), &[AlleleId::new(2)]);
        assert_eq!(idx.gene_of(AlleleId::new(1)), g12);
        assert!(idx.gene_id("nope").is_none());
    }
