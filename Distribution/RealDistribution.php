<?php

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

/**
 * Base interface for distributions on the reals.
 *
 * This class has been ported from Java: 
 *
 *    package org.apache.commons.math3.distribution.RealDistribution
 */
interface RealDistribution 
{
    /**
     * For a random variable {@code X} whose values are distributed according
     * to this distribution, this method returns {@code P(X = x)}. In other
     * words, this method represents the probability mass function (PMF)
     * for the distribution.
     *
     * @param x the point at which the PMF is evaluated
     * @return float
     */
    public function probability($x);

    /**
     * Returns the probability density function (PDF) of this distribution
     * evaluated at the specified point {@code x}. In general, the PDF is
     * the derivative of the {@link #cumulativeProbability(double) CDF}.
     * If the derivative does not exist at {@code x}, then an appropriate
     * replacement should be returned, e.g. {@code Double.POSITIVE_INFINITY},
     * {@code Double.NaN}, or  the limit inferior or limit superior of the
     * difference quotient.
     *
     * @param x the point at which the PDF is evaluated
     * @return float the value of the probability density function at point {@code x}
     */
    public function density($x);

    /**
     * For a random variable {@code X} whose values are distributed according
     * to this distribution, this method returns {@code P(X <= x)}. In other
     * words, this method represents the (cumulative) distribution function
     * (CDF) for this distribution.
     *
     * @param x the point at which the CDF is evaluated
     * @return float the probability that a random variable with this
     * distribution takes a value less than or equal to {@code x}
     */
    public function cumulativeProbability($x);

    /**
     * Computes the quantile function of this distribution. For a random
     * variable {@code X} distributed according to this distribution, the
     * returned value is
     * <ul>
     * <li><code>inf{x in R | P(X<=x) >= p}</code> for {@code 0 < p <= 1},</li>
     * <li><code>inf{x in R | P(X<=x) > 0}</code> for {@code p = 0}.</li>
     * </ul>
     *
     * @param p the cumulative probability
     * @return float the smallest {@code p}-quantile of this distribution
     * (largest 0-quantile for {@code p = 0})
     * @throws \OutOfRangeException if {@code p < 0} or {@code p > 1}
     */
    public function inverseCumulativeProbability($p);

    /**
     * Use this method to get the numerical value of the mean of this
     * distribution.
     *
     * @return the mean or {@code Double.NaN} if it is not defined
     */
    public function getNumericalMean();

    /**
     * Use this method to get the numerical value of the variance of this
     * distribution.
     *
     * @return float the variance (possibly {@code Double.POSITIVE_INFINITY} as
     * for certain cases in {@link TDistribution}) or {@code Double.NaN} if it
     * is not defined
     */
    public function getNumericalVariance();

    /**
     * Access the lower bound of the support. This method must return the same
     * value as {@code inverseCumulativeProbability(0)}. In other words, this
     * method must return
     * <p><code>inf {x in R | P(X <= x) > 0}</code>.</p>
     *
     * @return float lower bound of the support (might be
     * {@code Double.NEGATIVE_INFINITY})
     */
    public function getSupportLowerBound();

    /**
     * Access the upper bound of the support. This method must return the same
     * value as {@code inverseCumulativeProbability(1)}. In other words, this
     * method must return
     * <p><code>inf {x in R | P(X <= x) = 1}</code>.</p>
     *
     * @return float upper bound of the support (might be
     * {@code Double.POSITIVE_INFINITY})
     */
    public function getSupportUpperBound();

    /**
     * Use this method to get information about whether the support is connected,
     * i.e. whether all values between the lower and upper bound of the support
     * are included in the support.
     *
     * @return boolean whether the support is connected or not
     */
    public function  isSupportConnected();

    /**
     * Reseed the random generator used to generate samples.
     *
     * @param seed the new seed
     * @return void
     */
    public function reseedRandomGenerator($seed);

    /**
     * Generate a random value sampled from this distribution.
     *
     * @return float a random value.
     */
    public function sample();

    /**
     * Generate a random sample from the distribution.
     *
     * @param sampleSize the number of random values to generate
     * @return float[] an array representing the random sample
     * @throws org.apache.commons.math3.exception.NotStrictlyPositiveException
     * if {@code sampleSize} is not positive
     */
    public function sample($sampleSize);
}
