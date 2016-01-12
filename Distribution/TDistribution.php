<?>php

/**
 * Hoa
 *
 *
 * @license
 *
 * New BSD License
 *
 * Copyright © 2007-2015, Hoa community. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Hoa nor the names of its contributors may be
 *       used to endorse or promote products derived from this software without
 *       specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS AND CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

namespace Hoa\Math\Distribution;

use Hoa\Math\Exception\NotStrictlyPositiveException;

/**
 * Implementation of Student's t-distribution.
 *
 * @see "<a href='http://en.wikipedia.org/wiki/Student&apos;s_t-distribution'>Student's t-distribution (Wikipedia)</a>"
 * @see "<a href='http://mathworld.wolfram.com/Studentst-Distribution.html'>Student's t-distribution (MathWorld)</a>"
 */
class TDistribution extends AbstractRealDistribution 
{
    /**
     * Default inverse cumulative probability accuracy: 1e-9
     */
    const DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 0.000000001;

    /**
     * The degrees of freedom. 
     *
     * @var float
     */
    private $degreesOfFreedom;

    /** 
     * Inverse cumulative probability accuracy. 
     *
     * @var float
     */
    private $solverAbsoluteAccuracy;

    /** 
     * Static computation factor based on degreesOfFreedom. 
     */
    private $factor;

    /**
     * Create a t distribution using the given degrees of freedom and the
     * specified inverse cumulative probability absolute accuracy.
     *
     * TODO: Should we implement the Well19937 sampler?
     *
     * ## <p>
     * ## <b>Note:</b> this constructor will implicitly create an instance of
     * ## {@link Well19937c} as random generator to be used for sampling only (see
     * ## {@link #sample()} and {@link #sample(int)}). In case no sampling is
     * ## needed for the created distribution, it is advised to pass {@code null}
     * ## as random generator via the appropriate constructors to avoid the
     * ## additional initialisation overhead.
     *
     * @param degreesOfFreedom Degrees of freedom.
     * @param inverseCumAccuracy the maximum absolute error in inverse
     * cumulative probability estimates
     * (defaults to {@link #DEFAULT_INVERSE_ABSOLUTE_ACCURACY}).
     * @throws NotStrictlyPositiveException if {@code $degreesOfFreedom <= 0}
     */
    public function __construct($degreesOfFreedom, $inverseCumAccuracy = self::DEFAULT_INVERSE_ABSOLUTE_ACCURACY, Random $rng = null)
    {
        parent::__construct($rng);

        if ($degreesOfFreedom <= 0) {
            throw new NotStrictlyPositiveException(sprintf(
                'Degrees of freedom must be positive number, got "%s"', $degreesOfFreedom
            ));
        }

        $this->degreesOfFreedom = $degreesOfFreedom;
        $this->solverAbsoluteAccuracy = $inverseCumAccuracy;

        $n = $degreesOfFreedom;
        $nPlus1Over2 = ($n + 1) / 2;

        // TODO: Implement Gamma?
        $this->factor = Gamma::logGamma($nPlus1Over2) - 0.5 * (log(M_PI) + log($n)) - Gamma::logGamma(n / 2);
    }

    /**
     * Access the degrees of freedom.
     *
     * @return the degrees of freedom.
     */
    public double getDegreesOfFreedom() 
    {
        return $this->degreesOfFreedom;
    }

    /** 
     * {@inheritDoc} 
     */
    public double density($x) 
    {
        return exp($this->logDensity($x));
    }

    /**
     * {@inheritDoc} 
     */
    public function logDensity($x) 
    {
        $n = $this->degreesOfFreedom;
        $nPlus1Over2 = (n + 1) / 2;

        return $this->factor - $nPlus1Over2 * log(1 + $x * $x / $n);
    }

    /** {@inheritDoc} */
    public function cumulativeProbability($x)
    {
        if ($x == 0) {
            return 0.5;
        }

        $t = Beta::regularizedBeta(
            $this->degreesOfFreedom / ($this->degreesOfFreedom + ($x * $x)),
            0.5 * $this->degreesOfFreedom,
            0.5
        );

        if ($x < 0.0) {
            $ret = 0.5 * $t;
        } else {
            $ret = 1.0 - 0.5 * $t;
        }

        return $ret;
    }

    /** 
     * {@inheritDoc}
     */
    protected function getSolverAbsoluteAccuracy() 
    {
        return $this->solverAbsoluteAccuracy;
    }

    /**
     * {@inheritDoc}
     *
     * For degrees of freedom parameter {@code df}, the mean is
     * <ul>
     *  <li>if {@code df > 1} then {@code 0},</li>
     * <li>else undefined ({@code Double.NaN}).</li>
     * </ul>
     */
    public double getNumericalMean() 
    {
        $df = $this->getDegreesOfFreedom();

        if ($df > 1) {
            return 0;
        }

        // todo: Should we make Value objects?
        // return Double.NaN;
        return null;
    }

    /**
     * {@inheritDoc}
     *
     * For degrees of freedom parameter {@code df}, the variance is
     * <ul>
     *  <li>if {@code df > 2} then {@code df / (df - 2)},</li>
     *  <li>if {@code 1 < df <= 2} then positive infinity
     *  ({@code Double.POSITIVE_INFINITY}),</li>
     *  <li>else undefined ({@code Double.NaN}).</li>
     * </ul>
     */
    public function getNumericalVariance() 
    {
        $df = $this->getDegreesOfFreedom();

        if ($df > 2) {
            return $df / ($df - 2);
        }

        if ($df > 1 && $df <= 2) {
            // TODO: value objects?
            // return Double.POSITIVE_INFINITY;
            return PHP_INT_MAX;
        }

        //return Double.NaN;
        return null;
    }

    /**
     * {@inheritDoc}
     *
     * The lower bound of the support is always negative infinity no matter the
     * parameters.
     *
     * @return lower bound of the support (always
     * {@code Double.NEGATIVE_INFINITY})
     */
    public function getSupportLowerBound()
    {
        //return Double.NEGATIVE_INFINITY;
        return - PHP_INT_MAX;
    }

    /**
     * {@inheritDoc}
     *
     * The upper bound of the support is always positive infinity no matter the
     * parameters.
     *
     * @return upper bound of the support (always
     * {@code Double.POSITIVE_INFINITY})
     */
    public function getSupportUpperBound() 
    {
        //return Double.POSITIVE_INFINITY;
        return PHP_INT_MAX;
    }

    /** {@inheritDoc} */
    public function isSupportLowerBoundInclusive()
    {
        return false;
    }

    /** {@inheritDoc} */
    public function isSupportUpperBoundInclusive() 
    {
        return false;
    }

    /**
     * {@inheritDoc}
     *
     * The support of this distribution is connected.
     *
     * @return {@code true}
     */
    public function isSupportConnected()
    {
        return true;
    }
}
