import pytest
from scripts.insert_virus import Events

class TestEvent:

    @ pytest.mark.parametrize("args, kwargs, expected", [
        ((), {}, TypeError),
        (("", ""), {}, OSError),
        (("test", "test"), {}, OSError),
        (("test/refs/test_human.fa", "test/refs/test_virus.fa"), 
            {'fasta_extensions':"blah"}, AssertionError),
        ],
    )
    def test_init_bad_inputs(self, args, kwargs, expected):
        with pytest.raises(expected):
            Events(*args, **kwargs)
