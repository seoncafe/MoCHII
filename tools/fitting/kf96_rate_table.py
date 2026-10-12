"""Resolve published KF96 charge-transfer fits on a log-log rate grid.

The existing Fortran CXRTAB evaluator linearly interpolates log(rate) against
log(temperature). Adaptive refinement preserves the original KF96 fit over
its own published temperature interval; the end nodes define clamping.
"""

from bisect import bisect_right
import heapq
import math


def charge_transfer_rate_on_grid(temperature, temperatures, rates):
    """Numerically reproduce the Fortran CXRTAB log-log interpolation."""
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be positive and finite")
    if temperature <= temperatures[0]:
        return rates[0]
    if temperature >= temperatures[-1]:
        return rates[-1]
    i = bisect_right(temperatures, temperature) - 1
    fraction = (math.log(temperature) - math.log(temperatures[i])) / (
        math.log(temperatures[i + 1]) - math.log(temperatures[i]))
    return math.exp((1 - fraction) * math.log(rates[i]) +
                    fraction * math.log(rates[i + 1]))


def kf96_charge_transfer_table(reaction, max_nodes=1024):
    """Return `(temperatures, rates)` for one published positive KF96 fit.

    The worst error is sampled at the quarter, middle, and three-quarter
    positions in each log-temperature segment. A second independent 2001-point
    check must be performed by the caller after output decimal formatting.
    `max_nodes` is a finite stop condition and should not exceed the Fortran
    reader's private MAX_CXR_TAB setting.
    """
    if max_nodes < 2:
        raise ValueError("max_nodes must be at least two")
    start = math.log(reaction.min_temperature)
    end = math.log(reaction.max_temperature)
    if not (math.isfinite(start) and math.isfinite(end) and start < end):
        raise ValueError("invalid published temperature interval")

    def log_rate(x):
        value = reaction.rate(math.exp(x))
        if not math.isfinite(value) or value <= 0:
            raise ValueError("KF96 fit has a nonpositive rate")
        return math.log(value)

    def segment_error(left, right):
        first, last = log_rate(left), log_rate(right)
        return max(abs(math.expm1(
            first * (1 - f) + last * f -
            log_rate(left * (1 - f) + right * f)))
            for f in (0.25, 0.5, 0.75))

    # Refining to half the requested tolerance leaves room for printed-value
    # rounding and unsampled positions in a segment.
    target = 5e-6
    segments = [(-segment_error(start, end), start, end)]
    while -segments[0][0] > target:
        if len(segments) + 1 >= max_nodes:
            raise ValueError("MAX_CX_TAB is too small for KF96 fit accuracy")
        _, left, right = heapq.heappop(segments)
        middle = 0.5 * (left + right)
        heapq.heappush(segments,
                       (-segment_error(left, middle), left, middle))
        heapq.heappush(segments,
                       (-segment_error(middle, right), middle, right))
    nodes = sorted({start, end, *(left for _, left, _ in segments),
                     *(right for _, _, right in segments)})
    temperatures = [reaction.min_temperature if x == start else
                    reaction.max_temperature if x == end else math.exp(x)
                    for x in nodes]
    rates = [reaction.rate(t) for t in temperatures]
    return temperatures, rates
