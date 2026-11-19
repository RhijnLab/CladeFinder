import pytest
from cladefinder.phylo import filter_variable_sites


def test_withheld_unique_alleles_become_invariant():
    # Positions 1 and 4 vary among references, but position 3 only varies
    # when the withheld isolate is included; it must not be retained.
    refs = {'a': 'AACA', 'b': 'AGCA', 'c': 'AGCT'}
    query = 'ATGA'
    filtered, q, stats = filter_variable_sites(refs, query)
    assert stats == {'original_sites': 4, 'retained_sites': 2,
                     'removed_invariant_or_ambiguous_sites': 2}
    assert filtered == {'a': 'AA', 'b': 'GA', 'c': 'GT'}
    assert q == 'TA'
    assert all(len(v) == len(q) for v in filtered.values())


def test_invariant_sites_after_tip_removal():
    refs = {'a': 'AAAA', 'b': 'AACA', 'c': 'AACA'}
    filtered, q, stats = filter_variable_sites(refs, 'GGGG')
    assert stats['retained_sites'] == 1
    assert list(filtered.values()) == ['A', 'C', 'C']
    assert q == 'G'


def test_ambiguous_training_columns_excluded():
    refs = {'a': 'AN', 'b': 'CN', 'c': 'GN'}
    filtered, q, stats = filter_variable_sites(refs, 'TN')
    assert stats['retained_sites'] == 1
    assert q == 'T'


def test_empty_or_mismatched_rejected():
    with pytest.raises(ValueError, match='lengths disagree'):
        filter_variable_sites({'a': 'AC', 'b': 'AT'}, 'A')
    with pytest.raises(ValueError, match='No variable'):
        filter_variable_sites({'a': 'AAAA', 'b': 'AAAA'}, 'AAAA')
