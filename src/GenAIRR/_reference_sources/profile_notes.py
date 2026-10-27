"""Span-limiting curation notes of the validated OGRDB exports, read by hand (data only).

Every curation note of the exports listed in the interpretation profile's scope that limits an
allele's affirmed span: an affirmed column ("up to and including nucleotide 319" or "... base
319") or a trailing-dot marker for likely further bases.  Each row is

``(set id, release, label, stated column, curation text sha256, airr coding sha256,
airr_ex coding sha256, verbatim sentences)``

The stated column is the last affirmed column as the note states it (``None`` when the note
only marks likely further bases).  The sentences are copied exactly from the note and include
the provider's own confirmations.  The text and coding sequences are those of both record
exports of the allele.  Nothing here is parsed or interpreted at run time.
"""

NOTES = (
    ("9606.IGH_VDJ", 10, "IGHV1-2*06", 319,
     "50e0c3f5f0464a2ee71b91ab91aec77baefc65e70e823752f3e1c7aa35d7b893",
     "5cf10cbedfd994b1da08a1ba73b261238bc48fdbf815be33c077f78342acd34b",
     "5cf10cbedfd994b1da08a1ba73b261238bc48fdbf815be33c077f78342acd34b", (
      "The full length sequence was first reported to IgPdb in 2015 by Scheepers and colleagues "
      "(see J Immunol 194:4371-8) from genomic sequencing, and it was separately reported to IgPdb "
      "in 2016 by Kirik and colleagues (see Mol Immunol 87:12-22), as an inference from VDJ "
      "rearrangements.",
      "In line with IARC policy, the submitted sequence was recognized up to and including "
      "nucleotide 319.",
      "(April 12, 2019) A trailing \".\" in the sequence indicates IARC's opinion that the sequence "
      "is likely to contain additional 3' nucleotides for which there is insufficient evidence to "
      "make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV1-3*04", None,
     "d59b28df76fbb94538c34880a7e30325abae4ef08dca2315dbe97aef905958ec",
     "61dd36cc5967fd86156f9f617641a1118d4b6e7216e4b35e389e95376416a0d7",
     "61dd36cc5967fd86156f9f617641a1118d4b6e7216e4b35e389e95376416a0d7", (
      "(April 12, 2019) A trailing \".\" in the sequence indicates IARC's opinion that the sequence "
      "is likely to contain additional 3' nucleotides for which there is insufficient evidence to "
      "make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV1-58*03", 319,
     "7c7a05c310246eac29d0da6d799c747279be970f2c236a524cf618b235de786a",
     "5088ac8013c8108c0d931aab0f1d1b6eeeba7e5e3d6cb21d68d2bcc2a0c96107",
     "ffd59d5d195ebdff16a299b296047d87d0df6cdf46caa25ace045ffa117c20d6", (
      "In line with previous policy, the submitted sequence will be recognized up to and including "
      "nucleotide 319.",
      "A trailing \u201c.\u201d indicates IARC's opinion that the sequences likely to contain an "
      "additional/additional 3' nucleotide(s) for which there is insufficient evidence to make an "
      "affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV1-69*18", 319,
     "93cde878d8a694ea8ba2b47c1129445193e2d90ddf3a1a402d36c80e6f59f657",
     "9217ea33e1a4572f71a74817260f3338bd919aee9f839e7251dbc2d697335b4a",
     "2c3d58af25ec6844a36fd9c4d8615e392d0546ee99113152283bf815daea7e60", (
      "There was strong support for the sequence being affirmed up to and including nucleotide 319, "
      "in line with IARC policy, as a Level 1 sequence, with one additional 3\u2019 nucleotide being "
      "likely present in the sequence.",
      "This is indicated in IARC and OGRDB publications by one dot at the end of the sequence.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV1-69*20", 319,
     "55d2ae4e4b6f4bd49f335909f62f76bc575c9e8f0ee1e2ee1b110d51a6bea063",
     "69785bea9144e5661fe4d1f2bbb088a28b443b17656d1509b0c2648d98f931ae",
     "5e4cbe8e4106af7fd4d657b3418cf4833c4fdd4e33d82e8a2014d461b475c4d5", (
      "IARC affirms the sequence at Level 1 only up to and including base 319 in agreement with "
      "past practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A, at "
      "base position 320.",
      "A trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain an additional "
      "3\u2019 nucleotide for which there is insufficient evidence to make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV1-69*i04", 319,
     "5bc86edac3e0addeb169bfe387b1bc78cb21c771428f721c38d13d8329b39add",
     "dc9ff3bf01bf29781c9306efada5620d17a57151c85937a140c4459c91dadf24",
     "dc9ff3bf01bf29781c9306efada5620d17a57151c85937a140c4459c91dadf24", (
      "IARC affirms the sequence at level 1 up to and including base 319 in agreement with past "
      "practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A at "
      "base position 320.",
      "A trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain one "
      "additional 3\u2019-nucleotides for which there is insufficient evidence to make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV2-70*20", 319,
     "7780dad70880b5228ecbf581ed1366fa35afcc204c5c145184256c0b5d45056b",
     "26318b815c64b8a61668b0efdc1d9cf58b14cf11c1c44a36a60bec6d7e6f3864",
     "26318b815c64b8a61668b0efdc1d9cf58b14cf11c1c44a36a60bec6d7e6f3864", (
      "The sequence was affirmed up to and including nucleotide 319 even though alleles of "
      "IGHV2-70/IGHV2-70D in general extend up to and including base 322.",
      "It is however likely that three additional 3\u2019 nucleotides are present in the sequence, and "
      "this is indicated in IARC and OGRDB publications by three dots at the end of the sequence.",
      "This sequence was incorporated into the IMGT germline gene database as a consequence of the "
      "IMGT-NC Report #2020-2-1008 (October 8th, 2020).",
      "The sequence that is now recognized as IGHV2-70*20 includes the final three nucleotides of "
      "the sequence - TAC.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-13*06", 319,
     "b5f824e8bc62d4229a45aecab6ec876f8d87d1f8bd6cf2efd5ffeb7993e3f251",
     "552c0d9bb4f9c105560536de573ba212a2efd54834b76b6745d1f3e83911959b",
     "552c0d9bb4f9c105560536de573ba212a2efd54834b76b6745d1f3e83911959b", (
      "IARC affirms the sequence as IGHV3-13*i01 at Level 1 up to and including base 319.",
      "It is acknowledged that the allele most likely carries 1 additional base, typically A, at "
      "base positions 320.",
      "Trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain additional "
      "3\u2019-nucleotides for which there is insufficient evidence to make an affirmation.",
      "For use in a reference germline gene set, IARC recommends the use of the expected full "
      "length sequence.",
      "The missing terminal nucleotide of the sequence can be confirmed as \u2018A\u2019 by reference to "
      "recent genomic sequences of Rodriguez and colleagues that are reported in the VDJbase "
      "database from eight individuals.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-21*05", 319,
     "a7c45550f6ea2b2b4b971c633fb5a33847da8bcea92458b3f81767a57dfae1c6",
     "35cb8ca3a1b160541e7077d73b204ce4aeff6ad242afc84adda52e5f1313c3ff",
     "95560b1a0447d3a41093360f9dbaa89679e024319ae01531cf06d502bd8be1f6", (
      "In line with previous policy, the submitted sequence is recognized up to and including "
      "nucleotide 319.",
      "A trailing \u201c.\u201d in the affirmed sequence indicates IARC\u2019s opinion that the sequence is likely "
      "to contain an additional 3\u2019 nucleotide for which there is insufficient evidence to make an "
      "affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-21*07", 320,
     "6d7213d0aa90ee5e2260070c971eae6f346d8832d779a2a4166bf5b677eeb2f9",
     "276e4c5846cff333ecb1f65c5b4258e3843d9dba3bc88eba3f6fe804e1a9e414",
     "276e4c5846cff333ecb1f65c5b4258e3843d9dba3bc88eba3f6fe804e1a9e414", (
      "The sequence was tentatively affirmed as a Level 1 sequence, and the final 3\u2019 nucleotides "
      "will be considered at a later date, when terminal nucleotide plots become available, at "
      "which time the affirmed sequence will be noted in the IARC minutes.",
      "IARC affirms the sequence up to and including base 320 based on the convincing evidence of "
      "the sequence end as presented in the OGRDBplot.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-30*20", 319,
     "91a7dd0d8cc1128dc85c0e5a9e1bc3fd898ffa07563a0b6a858ba849786ebb80",
     "a73fbf9be6ae086780a2f3bae56d3ca5830cbb3006ae464e77f0bcb2e861cbcc",
     "acdfb027f8c473f816fce4e894fa8bc9da1bb68cabf7575462d13ee67bc2d722", (
      "IARC affirms the sequence at Level 1 up to and including base 319 in agreement with past "
      "practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A at "
      "base position 320.",
      "A trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain an additional "
      "3\u2019 nucleotide for which there is insufficient evidence to make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-30*i02", 319,
     "f6b548c962c679a6911411170524dd0909f39dcc845f8b6d03702f2cca851834",
     "f9baa49c6ef688109d6c75b8b4dda01747265a3193350cda92e57f03929fe201",
     "f1bfdb678b164515233d82c81c55e0be85a1c919befad0fce0fdbc1beb77034b", (
      "IARC now affirms the sequence at level 1 up to and including base 319.",
      "It is acknowledged that the allele most likely carries 1 additional base, typically A, at "
      "base positions 320.",
      "The trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain one "
      "additional 3\u2019-nucleotide for which there is insufficient evidence to make an affirmation.",
      "The GenBank record (OX384049) is publicly available and reports the sequence of IGHV3-30*04 "
      "C201T G317A, up to and including base 320.",
      "For use in a reference germline gene set, IARC recommends the use of the expected full "
      "length sequence.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-33*08", 319,
     "8e3d6c364d3045603597260dbc3c02cee5a45f9c441eca32bad1a09e30ae634b",
     "9be09b70d7d51dc52263e1c56d12fda9ee2e1f4137de40bc9431c7bd46245aae",
     "9be09b70d7d51dc52263e1c56d12fda9ee2e1f4137de40bc9431c7bd46245aae", (
      "IGHV3-33*01_G75C is identical to IGHV3-33*08 that entered into the IMGT-DB on March 5th, "
      "2021 based on its presence in GenBank entry accession number AC279998.",
      "IARC affirms, based on information in the transcriptome of this subject, the sequence as a "
      "level 1 sequence up to and including base 319 in agreement with past practice.",
      "A trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain an additional "
      "3\u2019 nucleotide for which there is insufficient evidence to make an affirmation.",
      "It is acknowledged that the allele most likely carry one additional base, typically A at "
      "base position 320.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-43D*04", 321,
     "5f8f77e25afda457f5e02d5a330f71d3b58a6bae159425ecbf97dc7e5ec06615",
     "501c18e7a0fdee1f9928c012a8dd2d9234c453f73d756393cd7e43bca2708a12",
     "501c18e7a0fdee1f9928c012a8dd2d9234c453f73d756393cd7e43bca2708a12", (
      "(April 12, 2019) In line with IARC policy, the submitted sequence was recognized up to and "
      "including nucleotide 321.",
      "A trailing \u201c.\u201d indicates IARC's opinion that the sequences likely to contain an "
      "additional/additional 3' nucleotide(s) for which there is insufficient evidence to make an "
      "affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-43D*05", 321,
     "017ae3a05e22f6cdacf4078252961cf03ec49aff2ea3bd611326d7da3dab72d7",
     "9c3f83fc9e56c3c09079bf167bc1069ca9c92f6c0d378eaea0fa9e38d17b740a",
     "9c3f83fc9e56c3c09079bf167bc1069ca9c92f6c0d378eaea0fa9e38d17b740a", (
      "IARC affirms the sequence at level 1 up to and including base 321 in agreement with past "
      "practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A at "
      "base position 322.",
      "A trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain one "
      "additional 3\u2019-nucleotides for which there is insufficient evidence to make an affirmation.",
      "The missing terminal nucleotide of the sequence can be confirmed as \u2018A\u2019 by reference to "
      "recent genomic sequences of Rodriguez and colleagues that are reported in the VDJbase "
      "database from three individuals.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-53*05", 319,
     "e00ab4d3de805f66151b6b117532e53f3f3e26bf36f5c195b6f6e08751cf5747",
     "426f84866803db7ac8c601e5f92e006b8cb040ab223b9a8d9a5a9553d4bb7f8f",
     "cdd15f3aaf689925d30185433e22bc08e77055a747e2fc880fdd77a42b274b72", (
      "In line with previous policy, the submitted sequence was recognized up to and including "
      "nucleotide 319.",
      "A trailing \u201c.\u201d in the affirmed sequence indicates IARC\u2019s opinion that the sequence is likely "
      "to contain an additional 3\u2019 nucleotide for which there is insufficient evidence to make an "
      "affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-7*04", 319,
     "9f75b757a5e3397429a3040e2e42f9822487b68fa5034699411e25873afbc4a6",
     "7ca93cfaf7876b9f6088fa7f97558e08ea251118b3f5d8aa832dcefa30899acb",
     "69f2aa2b18222aea786a2e05ed1acc181107ab5c5b8430c35ba1f7a8c3ce0c0f", (
      "(April 12, 2019) In line with IARC policy, the submitted sequence was recognized up to and "
      "including nucleotide 319.",
      "A trailing \u201c.\u201d indicates IARC's opinion that the sequences likely to contain an "
      "additional/additional 3' nucleotide(s) for which there is insufficient evidence to make an "
      "affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV3-9*04", 319,
     "9d7566b908d22466458fc9f62e56ffa5ca3ddce62aa68184c45d4dc504a70f70",
     "3157f3f7e0676547f628a5562dcdb0b21f81068a8fd1162475abca2762108fb7",
     "c5204d3cd655b88d93dcf7c8daf262a95b1b2fcf78bbbe9f6f9004f6ef26bc05", (
      "IARC affirms the sequence up to and including base 319 as a level 1 sequence.",
      "It is acknowledged that the allele most likely carries 3 additional bases, typically ATA at "
      "base position 320-322.",
      "The trailing \u201c...\u201d indicates IARC\u2019s opinion that the sequence is likely to contain three "
      "additional 3\u2019 nucleotides for which there is insufficient evidence to make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-30-4*08", 319,
     "e03cb9370703d010d9a537c7eb59808b6fd03860f2a4cb79554becd0de778d37",
     "dcd1f9d83aa38a9f400477d3bb70edd98e5a87f7b815c0911e4780a6841726b0",
     "dcd1f9d83aa38a9f400477d3bb70edd98e5a87f7b815c0911e4780a6841726b0", (
      "The sequence has previously been reported from genomic sequencing, and by inference from VDJ "
      "data, and is documented in IgPdb as IGHV4-30-4*p08.",
      "In line with IARC policy, the submitted sequence was recognized up to and including "
      "nucleotide 319.",
      "(April 12, 2019) A trailing \".\" in the sequence indicates IARC's opinion that the sequence "
      "is likely to contain additional 3' nucleotides for which there is insufficient evidence to "
      "make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-30-4*i02", 319,
     "dd7f2696a8ba6d41cbed53b5a9987d3ada4e7ac35e0e0b5e9d5d66a4ba4b6f61",
     "4dcb55103a63b1ea88f3a0ea253c633440d65957e93335fcc879f78e809e38f0",
     "bb027ed661b75be9f773caffb07d1b9375bf1243aa1f9278576a259dbfd5e99a", (
      "IARC affirms the sequence at level 1 up to and including base 319 in agreement with past "
      "practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A at "
      "base position 320.",
      "A trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain one "
      "additional 3\u2019-nucleotides for which there is insufficient evidence to make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-31*11", 319,
     "b98ee70dc1e8d08a9e6bf8fa2db7ea50fbc5297addea6b8778c2b3c177097fa7",
     "a02ecba432123f9c08204431b91ca6a96a0c589311245c41dd2aeae2a5ac8114",
     "2bfb8ceec5b4c7c391c20a96d6d2e58e7cfceaedef7e934f139171d44a391f0b", (
      "In line with previous policy, the submitted sequence will be recognized up to and including "
      "nucleotide 319.",
      "(April 12, 2019) A trailing \".\" in the sequence indicates IARC's opinion that the sequence "
      "is likely to contain additional 3' nucleotides for which there is insufficient evidence to "
      "make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-39*08", 319,
     "7e5ac2ddd7569d18fd6d40a01db634bf6919b8419c605c0081c59f43b1d49959",
     "ff0902e040115105389fce57c169eb3ad91aeef1b0248bce5621eab5180d0e76",
     "a8e8166f52b52a9a8d3a9a253b43d93ce487d6be8c17f5432f205c337402adae", (
      "IARC affirms the sequence as a level 1 sequence up to and including base 319 in agreement "
      "with past practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A at "
      "base position 320.",
      "A trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain an additional "
      "3\u2019 nucleotide for which there is insufficient evidence to make an affirmation.",
      "The data describing the final bases of this inferred allele strongly supports that this base "
      "(319) is a C.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-39*09", 319,
     "0840331087c03e947083fcdde8ac5790a7ab67019ac43c1fde2070c3d4066fd1",
     "6764e2099114bb7470897b80ad7988ce23d4fc43712b20e7d11bd360a7fcde5f",
     "d9f2078bf716890a37d566fd90f757b7f7a630e193bc72339ae9ee770070ad03", (
      "The sequence was tentatively affirmed as a Level 1 sequence, and the final 3\u2019 nucleotides "
      "will be considered at a later date, when terminal nucleotide plots become available, at "
      "which time the affirmed sequence will be noted in the IARC minutes.",
      "IARC infers the sequence at level 1 up to and including base 319 in agreement with past "
      "practice.",
      "It is acknowledged that the allele most likely carry one additional base, typically A at "
      "base position 320.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-39*i03", 319,
     "6da2215223bd8b0de40e23c1b143142231e1180dcdcc35f3971cdec6e9b13456",
     "8193ac85d2ea62dc26418e14983dac755236190793f339885735e8a029f5b96c",
     "32fac52413e21f43e3e8a01145f7e1456389deca192562e2d797d4bcef4dff92", (
      "IARC affirms the sequence at level 1 up to and including base 319 in agreement with past "
      "practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A at "
      "base position 320.",
      "A trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain one "
      "additional 3\u2019-nucleotides for which there is insufficient evidence to make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-4*03", 319,
     "3a03d59f92f012e9b44aa33f70270638482bc8741120096510aa29852f8346ba",
     "5e19fb52419d2143c1694a560efaa286d96e301a1ec6bd9d55f719542714a795",
     "5e19fb52419d2143c1694a560efaa286d96e301a1ec6bd9d55f719542714a795", (
      "(April 12, 2019) In line with IARC policy, the submitted sequence was recognized up to and "
      "including nucleotide 319.",
      "A trailing \u201c.\u201d indicates IARC's opinion that the sequences likely to contain an "
      "additional/additional 3' nucleotide(s) for which there is insufficient evidence to make an "
      "affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-4*09", 318,
     "fd42b1c27fbb16048ffd75107dd0c67fb09bbd996355b4b5c4986d4b9b4e30ba",
     "33ff8619becba3b6a3e1b33a85d22390e41556eccdb757f5c2dc8b6f6ade7838",
     "f66abadca0c08bd80b47a6e6b9c4a24f9d2387eb73b3841c685a38f37748e31b", (
      "The sequence, up to and including nucleotide 318, was affirmed as a Level 1 sequence.",
      "Uncertainty regarding nucleotides 319 and 320 is indicated in IARC and OGRDB publications by "
      "two dots at the end of the affirmed sequence.",
      "The committee believes it is most likely that the 319th nucleotide is \u2018C\u2019.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-4*10", 319,
     "c0a259260ef44158334003951081dde2562ac7ee7e0659853cbdadeb55a733c2",
     "e8206d272eeb66b96a1a92cec2201c353f29597ebaa996cdc21b3c34870321d6",
     "dd8abdeca7bd75fc8c655a25eceeb0c35c5bbc6e6d0af8b888a751c5cbaa0010", (
      "IARC now affirms, based on the extensive validation, the sequence at level 1 up to and "
      "including base 319 in agreement with past practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A at "
      "base position 320.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-59*11", 319,
     "e660a3e581b71a1555a7b307034b40bc054e536ea4c9bbc2a1bb0ab1e899ec83",
     "b369293708e2e913bc2595e50375c88b25b3ed4a817e86bed81db79da79ccf20",
     "9529199174ba6576e76d61bbef262f86ab7a892e37e535c8452d5aa05eb7f7e6", (
      "The sequence, up to and including nucleotide 319, was affirmed as a Level 1 sequence, but "
      "was subsequently identified as the recently reported allele IGHV4-59*11.",
      "The sequence submitted to Genbank lacked the 320th nucleotide that is noted in the IMGT "
      "sequence.",
      "In IARC and OGRDB communications, a dot at the 3' end of the sequence shows that the IARC "
      "believes there is a nucleotide missing from the sequence.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-59*12", 319,
     "9391d61f6dc5694bd6ec08c16996d029856682f43279ac394b4832aa7f5b6048",
     "3aab7400d8caf0faf40198533502a512f8be521208fe7821b89296ada4c3b64d",
     "7eaf207f5741fedb91837c4fbcbdf73478433705b4083a91fa4ed2229a79f5ae", (
      "In line with previous policy, the submitted sequence will be recognized up to and including "
      "nucleotide 319.",
      "(April 12, 2019) A trailing \".\" in the sequence indicates IARC's opinion that the sequence "
      "is likely to contain additional 3' nucleotides for which there is insufficient evidence to "
      "make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-59*i03", 319,
     "57bfe95945a47a5b012b528f50c5ca57d56ba67ea9119bc188be92dd1aaea46b",
     "fc84af7ff6ee4a60309257c734c59d801c4bfdcf84cc82b1b7fe0a6806d87a18",
     "2900fbf031e170d3cce9202715e43f3b8628e885549a3e5cb001ce545fe33e58", (
      "IARC affirms the sequence at level 1 up to and including base 319 in agreement with past "
      "practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A at "
      "base position 320.",
      "A trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain one "
      "additional 3\u2019-nucleotides for which there is insufficient evidence to make an affirmation.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-61*10", 319,
     "6cac559941dad2dd9a0ce82d4c86a45170089a92eabbfb833e7f7d6b36e87092",
     "05a7b1e32f7a401c539cb87e15e1304b67cfa209c550c8c661c65b0cbd14d32b",
     "f131b64f8b627a27bd28573b0fd4f4446488062763611dce87416cf0fc21ea35", (
      "There was strong support for the sequence being affirmed up to and including nucleotide 319, "
      "in line with IARC policy, as a Level 1 sequence, with one additional 3\u2019 nucleotide being "
      "likely present in the sequence.",
      "This is indicated in IARC and OGRDB publications by one dot at the end of the sequence.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-61*11", 319,
     "8605e82d096851103a96d2600f91d5b613d079dedde59a9aef41ba4579e5f9ae",
     "25c8a3957692cf9b06b5dabff1aa70af11088fbb1f01ac33e29d56f6071fbdad",
     "c5399d6d6675cdd2c3e4386466a503b916309922a0c70810f495603f4751ca87", (
      "IARC, taking all of the collective evidence under consideration, affirms the sequence at "
      "level 1 up to and including base 319 in agreement with past practice.",
      "It is acknowledged that the allele most likely carries one additional base, typically A at "
      "base position 320.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV4-61*i04", 319,
     "3978b604b9ccd7b69e28d1bb4bfafa4a9e9ff7a0bf9b1b93c2d13326a5d18fb5",
     "355ba4225ca8fa74bcc9426636d7de7d9d45b123e250a067aa9294c1593675fd",
     "246012a8732ef6d583b077219bc900e7f232d61aa17993254fea558288e04d91", (
      "IARC affirms the sequence as IGHV4-61*i04 at Level 1 up to and including base 319.",
      "It is acknowledged that the allele most likely carries 1 additional base, typically A, at "
      "base positions 320.",
      "Trailing \u201c.\u201d indicates IARC\u2019s opinion that the sequence is likely to contain additional "
      "3\u2019-nucleotides for which there is insufficient evidence to make an affirmation.",
      "For use in a reference germline gene set, IARC recommends the use of the expected full "
      "length sequence.",
     )),
    ("9606.IGH_VDJ", 10, "IGHV5-51*06", 318,
     "583d3ea46e1f3f5e1e42d743ad4aa206e0bd53ad964a23fc905a3096c555a093",
     "f78037bfb74df2768a02f6f9517bcb65df39a39fe33140f1528a589819299929",
     "edb650ca13f27794bae615618eca6875f612c24286760a88b30677e1206fd000", (
      "The sequence is affirmed only up to and including base 318, the final base of the inference "
      "and its GenBank record, and of the corresponding closest allele (IGHV5-51*03) in the IMGT "
      "database.",
      "(April 12, 2019) A trailing \".\" indicates IARC's opinion that the sequence is likely to "
      "contain additional 3' nucleotides for which there is insufficient evidence to make an "
      "affirmation.",
     )),
    ("9606.IGKappa_VJ", 5, "IGKV1D-8*04", 332,
     "72cc9c80f0900a04d3eb36b0db519930125ccce3f1b9edbef13848b64c5cc420",
     "cb6e144ac93621ded35c6c2f8f95b12656cd03add1bcce3e9e18f4cc4f35f182",
     "9b5a400ceff1faa08693b48a3f63d8ca602a33c410cb7cacc89b406b7ffdcdbf", (
      "Plots of the final 3\u2019 nucleotides showed a high level of variability, making it impossible "
      "to determine the final three nucleotides with certainty.",
      "The sequence, up to and including nucleotide 332, was affirmed as a Level 1 sequence.",
      "Uncertainty regarding nucleotides 333-335 will be indicated in IARC and OGRDB publications "
      "by three dots at the end of the affirmed sequence.",
     )),
    ("9606.IGLambda_VJ", 4, "IGLV10-54*04", 339,
     "12408e2acae180953310dd6010927ba2abd467b8d556556d97f9c2e1504a5f70",
     "054d5b997ca16c9313b56e814b1604ce669dd2d5d42b8d93922a7686d1b8eec1",
     "03800d1a35cc28cd4586817b201cb54f09dbefe54b17dd15e73ca059675d7c68", (
      "Plots of the final 3\u2019 nucleotides showed a high level of variability, making it impossible "
      "to determine the final three nucleotides with certainty.",
      "The sequence, up to and including nucleotide 339, was affirmed as a Level 1 sequence.",
      "Uncertainty regarding nucleotides 340-341 will be indicated in IARC and OGRDB publications "
      "by two dots at the end of the affirmed sequence.",
     )),
    ("9606.IGLambda_VJ", 4, "IGLV2-14*03", 337,
     "3c06429abcbabf8c891cb8c95db7e2970060b8f4e94a60afe38d61220e20dafe",
     "5ca70268a23d9f1ba3e8381a3d9d1cd95743f6f6e3a7a7136057305e8d0db9d9",
     "a09488a27413bd19d5366eedd5c07dfe124e7ec952af381551427d1e1690d71e", (
      "There was support for the sequence being affirmed up to and including nucleotide 337, in "
      "line with IARC policy not to infer bases substantially affected by trimming during the "
      "rearrangement process, as a Level 1 sequence, with two additional 3\u2019 nucleotide being likely "
      "present in the sequence.",
      "This is indicated in IARC and OGRDB publications by two dots at the end of the sequence.",
      "A genomic sequence (MK587524) has following the inference been cloned from the genotype of "
      "the donor that confirms the inference.",
      "It also confirms that bases 338-339 that cannot be inferred with confidence of the gene are "
      "T and C, respectively.",
     )),
    ("9606.IGLambda_VJ", 4, "IGLV2-14*04", 337,
     "e7e7ea1f8b424ed5e93fa6d80279dd3d4b16fbdc3dad2a67880969d72baa2289",
     "9008dba762d2100cfb826c99c2cd5f173081bed025b4807a3d1d4e7c3384f061",
     "9008dba762d2100cfb826c99c2cd5f173081bed025b4807a3d1d4e7c3384f061", (
      "There was support for the sequence being affirmed up to and including nucleotide 337, in "
      "line with IARC policy not to infer bases substantially affected by trimming dring the "
      "rearrangement process, as a Level 1 sequence, with two additional 3\u2019 nucleotide being likely "
      "present in the sequence.",
     )),
    ("9606.IGLambda_VJ", 4, "IGLV3-21*04", 339,
     "d72b8ac22784b8260deb5bcccb4b49de0f34e7a4424b727aac6e69b0c3b0044e",
     "d0b983b1223e55be23ae8339d1e5060814900433fa90fd9a8d1e7cb4f3c8d285",
     "c1207f8ac1954e1610e6b6dbf0b7bc38301594c032dfc32ce293dc18f60c2b92", (
      "Plots of the final 3\u2019 nucleotides were considered, but the variability seen made it "
      "impossible to consider the final two nucleotides of the sequence.",
      "The sequence, up to and including nucleotide 339, was affirmed as the Level 1 sequence, "
      "IGLV3-21*i01.",
      "Uncertainty regarding nucleotides 340-341 will be indicated in IARC and OGRDB publications "
      "by two dots at the end of the affirmed sequence.",
     )),
    ("9606.IGLambda_VJ", 4, "IGLV3-25*03", 339,
     "73959a3aa995cbb0cb6f807584e709cf441f3745f0a06d8f2252196e5b80ad44",
     "47c8e1d37a1beed3957c0d059baf17113bd146f473d40ba4cd46b270dba43213",
     "59b875956dd72a662d02f278efb5a0f9537252de8ee84f09f6cb57d76e671d44", (
      "Plots of the final 3\u2019 nucleotides were considered, but the variability seen made it "
      "impossible to consider the final two nucleotides of the sequence.",
      "The sequence, up to and including nucleotide 339, was affirmed as the Level 1 sequence, "
      "IGLV3-25*i01.",
      "It appears to represent the full length version of the previously truncated IGLV3-25*03 "
      "sequence.",
      "It was subsequently noted that since the submission of this sequence, an identical genomic "
      "sequence has been accepted by IMGT as the full length sequence IGLV3-25*03.",
     )),
    ("9606.IGLambda_VJ", 4, "IGLV6-57*04", 333,
     "ba5d8ad8612f577e490f97bf07f874f3c4414d93454ec2bc602beb08259e57a3",
     "b7961ecc7ec1e4cd2707b66aaeb34e6bc24cb62e25f6f14245eb0ad3d83fe7eb",
     "9e9dcda340a220c8a7e4a15d71b9640b520ca1c6d2350fb35e5f989f6c6b886d", (
      "Plots of the final 3\u2019 nucleotides were considered, but the variability seen made it "
      "impossible to be certain of the final two nucleotides of the germline sequence.",
      "The sequence, up to and including nucleotide 333, was affirmed as the Level 1 sequence, "
      "IGLV6-57*i01.",
      "Uncertainty regarding nucleotides 334-335 will be indicated in IARC and OGRDB publications "
      "by two dots at the end of the affirmed sequence.",
     )),
)
