import numpy as np
import tensorflow as tf

from pytest import approx

from hmsc.utils.import_utils import calculate_GPP


SEED = 42


def input_values(rng):
    na = 2
    n1 = 4
    n2 = 3

    alpha = rng.random(na)
    d12 = rng.random(n1 * n2).reshape(n1, n2)
    d22 = rng.random(n2 * n2).reshape(n2, n2)
    d22 = 0.5 * (d22 + d22.T)
    np.fill_diagonal(d22, 0)
    return d12, d22, alpha


def reference_values():
    idD = \
[[ 4.74732221,  1.21841991, 16.81039263,  1.27917015],
 [ 2.89727739,  1.04729426,  3.83387891,  1.06836454]]
    iDW12 = \
[[[ 1.56552007,  1.92810424,  4.20341649],
  [ 0.3454143 ,  0.45571594,  0.44127379],
  [14.2458611 ,  9.3939915 , 10.41141587],
  [ 0.38626669,  0.55671448,  0.44182194]],
 [[ 0.40960001,  0.59143341,  2.33774095],
  [ 0.11340611,  0.18487531,  0.17466852],
  [ 2.86327373,  1.37390469,  1.64707451],
  [ 0.12931064,  0.24636518,  0.16388959]]]
    F = \
[[[13.80338756,  9.72259845, 10.89601792],
  [ 9.72259845,  7.44538439,  8.4114639 ],
  [10.89601792,  8.4114639 , 11.48249436]],
 [[ 3.22423038,  1.87731239,  1.82347872],
  [ 1.87731239,  1.70253007,  1.46121186],
  [ 1.82347872,  1.46121186,  3.64813777]]]
    iF = \
[[[ 0.90649386, -1.22934882,  0.04036143],
  [-1.22934882,  2.44625374, -0.62543626],
  [ 0.04036143, -0.62543626,  0.50695045]],
 [[ 0.88076752, -0.90416607, -0.07808988],
  [-0.90416607,  1.82323176, -0.27833386],
  [-0.07808988, -0.27833386,  0.4246276 ]]]
    detD = \
[-0.54608872, -0.15196945]
    return idD, iDW12, F, iF, detD


def test_calculate_GPP():
    rng = np.random.default_rng(seed=SEED)
    tf.keras.utils.set_random_seed(SEED)
    tf.config.experimental.enable_op_determinism()

    d12, d22, alpha = input_values(rng)

    values = calculate_GPP(d12, d22, alpha)
    values = list(map(lambda a: a.numpy(), values))
    names = ['idD', 'iDW12', 'F', 'iF', 'detD']
    assert len(names) == len(values)

    # Print values
    print()
    for name, array in zip(names, values):
        print(f'    {name} = \\')
        print(np.array2string(array, separator=', ', max_line_width=200))

    # Test against reference
    ref_values = list(map(np.asarray, reference_values()))
    assert len(ref_values) == len(values)
    for name, val, ref in zip(names, values, ref_values):
        assert val == approx(ref), name
