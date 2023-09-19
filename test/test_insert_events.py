import pytest
from conftest import hum_fa, vir_fa
from scripts.insert_virus import Events


@pytest.fixture
def events():
    return Events(hum_fa, vir_fa)

@pytest.fixture
def prob_dict():
    probs = {
        'p_whole': 0.5,
        'p_rearrange': 0.1,
        'p_delete': 0.1,
        'lambda_split': 5,
        'p_overlap': 0.1,
        'p_gap': 0.1,
        'lambda_junction': 5,
        'p_host_del': 0.5,
        'lambda_host_del': 10
    }
    return probs

class TestEvents:

    # should error because something wrong with args
    @ pytest.mark.parametrize("args, kwargs, expected", [
        ((), {}, TypeError),
        (("", ""), {}, OSError),
        (("test", "test"), {}, OSError),
        ((hum_fa, vir_fa), 
            {'seed':[]}, AssertionError),
        ((hum_fa, vir_fa), 
            {'seed':-1}, AssertionError),
        ((hum_fa, vir_fa), 
            {'seed':-1}, AssertionError),
        ((hum_fa, vir_fa), 
            {'min_len':-1}, AssertionError),
        ((hum_fa, vir_fa), 
            {'min_len':0}, AssertionError),
        ((hum_fa, vir_fa), 
            {'max_len':-1}, AssertionError),
        ((hum_fa, vir_fa), 
            {'max_len':0}, AssertionError),
        ((hum_fa, vir_fa), 
            {'max_len':10000000000000}, AssertionError),
        ((hum_fa, vir_fa), 
            {'max_len':10, 'min_len':20}, ValueError),
        ],
    )
    def test_init_bad_inputs(self, args, kwargs, expected):
        with pytest.raises(expected):
            Events(*args, **kwargs)

    # test of initialization
    @ pytest.mark.parametrize("args, kwargs", [
        ((hum_fa, vir_fa),
            {}),
        ])
    def test_init(self, args, kwargs):
        events = Events(*args, **kwargs)
        assert events is not None


    # test add_integrations should error
    @ pytest.mark.parametrize("num, expected", [
        (tuple(), TypeError),  # no probs or num
        ((0, ), AssertionError),  # zero 
        ((-1, ), AssertionError),  # negative num
        (("blah", ), AssertionError),  # string
        ((10000000000, ), ValueError),  # won't fit
    ])
    def test_add_integrations_bad_num(self, events, prob_dict, num, expected):
        with pytest.raises(expected):
            events.add_integrations(prob_dict *num)


    # test of add_integrations
    def test_add_integrations(self, events, prob_dict):
        events.add_integrations(prob_dict, 4)
        assert len(events.ints) == 4

