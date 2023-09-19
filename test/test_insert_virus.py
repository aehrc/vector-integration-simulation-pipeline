import pytest
from scripts.insert_virus import Events

class TestEvent:

    @ pytest.mark.parametrize("args, kwargs, expected", [
        ((), {}, TypeError),
        (("", ""), {}, OSError)]
    )
    def test_init_bad_inputs(self, args, kwargs, expected):
        with pytest.raises(expected):
            Events(*args, **kwargs)
