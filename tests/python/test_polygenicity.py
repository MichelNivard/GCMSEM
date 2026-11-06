import copy
import json
import numpy as np
import pytest
from gcmsem.polygenicity import summarize_polygenicity, json_safe


def fixture():
    # Two equally sized bins, identity LD. A point-normal model with pi=.1
    # has E(beta^2)=.1 and k4=3*pi*(1-pi)=.27, so M_e=M*pi=10.
    theta = np.array([[0., .1, 0., .27], [0., .1, 0., .27]])
    reps = np.repeat(theta[None, :, :], 4, axis=0)
    reps[:, :, 3] *= np.array([.8, .95, 1.05, 1.2])[:, None]
    blocks = [dict(chrom='1', start_row=i*25, stop_row=(i+1)*25) for i in range(4)]
    fit = dict(conditional=theta, conditional_jackknife=reps,
        metadata=dict(traits=['A'], jackknife_blocks=blocks, jackknife_fractions=[.25]*4))
    bins = np.tile([0, 1], 50)
    scores = dict(counts=np.array([50, 50]), af=np.where(bins == 0, .1, .4),
        L2=np.eye(2)[bins], bp=np.arange(1, 101), metadata=dict(maf_edges=[0, .25, .5]))
    return fit, scores


def test_point_normal_identity_and_covariance():
    fit, scores = fixture(); out = summarize_polygenicity(fit, scores)
    assert out['summary']['effective_count'] == pytest.approx(10)
    assert out['summary']['no_ld_count'] == pytest.approx(10)
    assert out['stages'][8]['estimate'] == pytest.approx(27)
    # Independent textbook equal-delete-group SE on the nonlinear ratio.
    r = np.array([row['effective_count'] for row in out['jackknife']])
    assert out['summary']['se'] == pytest.approx(np.sqrt(.75*np.sum((r-r.mean())**2)))
    assert sum(row['variance_fraction'] for row in out['jackknife']) == pytest.approx(1)
    # Perfectly correlated conditional perturbations require cross-bin terms.
    Q = np.asarray(out['conditional_covariance'])
    assert Q[1, 3] > 0


def test_scale_invariance_and_ld_effect():
    fit, scores = fixture(); base = summarize_polygenicity(fit, scores)
    altered = copy.deepcopy(fit)
    altered['conditional'] *= np.array([7., 49., 343., 2401.])
    altered['conditional_jackknife'] *= np.array([7., 49., 343., 2401.])
    out = summarize_polygenicity(altered, scores)
    for k in ('effective_count', 'se', 'delta_se', 'inverse_count', 'no_ld_count'):
        assert out['summary'][k] == pytest.approx(base['summary'][k])
    scores['L2'] *= 2
    out = summarize_polygenicity(fit, scores)
    assert out['summary']['effective_count'] == pytest.approx(300/33)
    assert out['summary']['no_ld_count'] == pytest.approx(10)


def test_singular_ratio_is_not_dropped_and_signed_values_survive():
    fit, scores = fixture()
    fit['conditional_jackknife'][0, :, 1] = 1.
    fit['conditional_jackknife'][0, :, 3] = -3.
    out = summarize_polygenicity(fit, scores)
    assert len(out['jackknife']) == 4
    assert 'nonpositive_denominator_in_deletions' in out['warnings']
    assert 'undefined_ratio_or_uncertainty' in out['warnings']
    assert not np.isfinite(out['jackknife'][0]['effective_count'])
    assert np.isnan(out['summary']['se'])
    assert json_safe(out)['summary']['se'] is None
    fit['conditional'][:, 3] = -.3
    out = summarize_polygenicity(fit, scores)
    assert out['summary']['effective_count'] < 0
    assert 'negative_fourth_cumulant_gaussian_mixture_incompatible' in out['warnings']
    json.dumps(json_safe(out), allow_nan=False)


def test_sparse_stratum_rejected():
    fit, scores = fixture()
    scores['af'][:25] = .1; scores['af'][25:] = .4
    with pytest.raises(ValueError, match='empties'):
        summarize_polygenicity(fit, scores)


def test_cli_protects_polygenicity_file_before_loading_inputs(tmp_path):
    from gcmsem.cli import parser, estimate_cli
    prefix = tmp_path / 'analysis'
    output = tmp_path / 'analysis.polygenicity.json'
    output.write_text('previous result\n')
    args = parser().parse_args(['polygenicity', '--sumstats', 'missing.tsv',
        '--scores', 'missing.npz', '--out', str(prefix)])
    with pytest.raises(FileExistsError, match='polygenicity output already exists'):
        estimate_cli(args)
    assert output.read_text() == 'previous result\n'
