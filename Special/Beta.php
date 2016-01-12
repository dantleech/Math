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

namespace Hoa\Math\Special;

/**
 * <p>
 * This is a utility class that provides computation methods related to the
 * Beta family of functions.
 * </p>
 * <p>
 * Implementation of {@link #logBeta(double, double)} is based on the
 * algorithms described in
 * <ul>
 * <li><a href="http://dx.doi.org/10.1145/22721.23109">Didonato and Morris
 *     (1986)</a>, <em>Computation of the Incomplete Gamma Function Ratios
 *     and their Inverse</em>, TOMS 12(4), 377-393,</li>
 * <li><a href="http://dx.doi.org/10.1145/131766.131776">Didonato and Morris
 *     (1992)</a>, <em>Algorithm 708: Significant Digit Computation of the
 *     Incomplete Beta Function Ratios</em>, TOMS 18(3), 360-373,</li>
 * </ul>
 * and implemented in the
 * <a href="http://www.dtic.mil/docs/citations/ADA476840">NSWC Library of Mathematical Functions</a>,
 * available
 * <a href="http://www.ualberta.ca/CNS/RESEARCH/Software/NumericalNSWC/site.html">here</a>.
 * This library is "approved for public release", and the
 * <a href="http://www.dtic.mil/dtic/pdf/announcements/CopyrightGuidance.pdf">Copyright guidance</a>
 * indicates that unless otherwise stated in the code, all FORTRAN functions in
 * this library are license free. Since no such notice appears in the code these
 * functions can safely be ported to Commons-Math.
 * </p>
 */
class Beta
{
    /** 
     * Maximum allowed numerical error. 1E-14
     */
    const DEFAULT_EPSILON = 1E-14;

    /** The constant value of ½log 2π. */
    const HALF_LOG_TWO_PI = .9189385332046727;

    /**
     * <p>
     * The coefficients of the series expansion of the Δ function. This function
     * is defined as follows
     * </p>
     * <center>Δ(x) = log Γ(x) - (x - 0.5) log a + a - 0.5 log 2π,</center>
     * <p>
     * see equation (23) in Didonato and Morris (1992). The series expansion,
     * which applies for x ≥ 10, reads
     * </p>
     * <pre>
     *                 14
     *                ====
     *             1  \                2 n
     *     Δ(x) = ---  >    d  (10 / x)
     *             x  /      n
     *                ====
     *                n = 0
     * <pre>
     */
    private $DELTA = [
        .833333333333333333333333333333E-01,
        -.277777777777777777777777752282E-04,
        .793650793650793650791732130419E-07,
        -.595238095238095232389839236182E-09,
        .841750841750832853294451671990E-11,
        -.191752691751854612334149171243E-12,
        .641025640510325475730918472625E-14,
        -.295506514125338232839867823991E-15,
        .179643716359402238723287696452E-16,
        -.139228964661627791231203060395E-17,
        .133802855014020915603275339093E-18,
        -.154246009867966094273710216533E-19,
        .197701992980957427278370133333E-20,
        -.234065664793997056856992426667E-21,
        .171348014966398575409015466667E-22
    ];

    /**
     * Default constructor.  Prohibit instantiation.
     */
    final private function __construct()
    {
    }

    /**
     * Returns the regularized beta function I(x, a, b).
     *
     * The implementation of this method is based on:
     * <ul>
     * <li>
     * <a href="http://mathworld.wolfram.com/RegularizedBetaFunction.html">
     * Regularized Beta Function</a>.</li>
     * <li>
     * <a href="http://functions.wolfram.com/06.21.10.0001.01">
     * Regularized Beta Function</a>.</li>
     * </ul>
     *
     * @param x the value.
     * @param a Parameter {@code a}.
     * @param b Parameter {@code b}.
     * @param epsilon When the absolute value of the nth item in the
     * series is less than epsilon the approximation ceases to calculate
     * further elements in the series.
     * @param maxIterations Maximum number of "iterations" to complete.
     * @return the regularized beta function I(x, a, b)
     * @throws org.apache.commons.math3.exception.MaxCountExceededException
     * if the algorithm fails to converge.
     */
    public static function regularizedBeta($x, $a, $b, $epsilon = self::DEFAULT_EPSILON, $maxIterations = PHP_INT_MAX)
    {
        $ret = null;

        if (NAN === $x ||
            NAN === $a ||
            NAN === $b ||
            $x < 0 ||
            $x > 1 ||
            $a <= 0 ||
            $b <= 0) {
            return NAN;
        }
       
        if ($x > ($a + 1) / (2 + $b + $a) &&
                   1 - $x <= ($b + 1) / (2 + $b + $a)) {
            return 1 - $this->regularizedBeta(1 - $x, $b, $a, $epsilon, $maxIterations);
        }


        $fraction = new ContinuedFraction(
            function () {
                return 1.0;
            },
            function ($n, $x) {
                $ret;
                $m;
                if ($n % 2 == 0) { // even
                    $m = $n / 2.0;
                    $ret = ($m * ($b - $m) * $x) /
                        (($a + (2 * $m) - 1) * ($a + (2 * $m)));
                } else {
                    $m = ($n - 1.0) / 2.0;
                    $ret = -(($a + $m) * ($a + $b + $m) * $x) /
                            (($a + (2 * $m)) * ($a + (2 * $m) + 1.0));
                }
                return $ret;
            }
        );

        $ret = exp(($a * log($x)) + ($b * log1p(-$x)) -
            log($a) - $this->logBeta($a, $b)) *
            1.0 / $fraction->evaluate($x, $epsilon, $maxIterations);

        return $ret;
    }


    /**
     * Returns the value of log Γ(a + b) for 1 ≤ a, b ≤ 2. Based on the
     * <em>NSWC Library of Mathematics Subroutines</em> double precision
     * implementation, {@code DGSMLN}. In {@code BetaTest.testLogGammaSum()},
     * this private method is accessed through reflection.
     *
     * @param float $a First argument.
     * @param float $b Second argument.
     * @return double The value of {@code log(Gamma(a + b))}.
     * @throws OutOfRangeException if {@code a} or {@code b} is lower than
     * {@code 1.0} or greater than {@code 2.0}.
     */
    private static function logGammaSum($a, $b)
    {
        if (($a < 1.0) || ($a > 2.0)) {
            throw new \OutOfRangeException($a, 1.0, 2.0);
        }
        if (($b < 1.0) || ($b > 2.0)) {
            throw new \OutOfRangeException($b, 1.0, 2.0);
        }

        $x = ($a - 1.0) + ($b - 1.0);
        if ($x <= 0.5) {
            return Gamma::logGamma1p(1.0 + $x);
        } else if (x <= 1.5) {
            return Gamma::logGamma1p($x) + log1p($x);
        } else {
            return Gamma::logGamma1p($x - 1.0) + log($x * (1.0 + $x));
        }
    }

    /**
     * Returns the value of log[Γ(b) / Γ(a + b)] for a ≥ 0 and b ≥ 10. Based on
     * the <em>NSWC Library of Mathematics Subroutines</em> double precision
     * implementation, {@code DLGDIV}. In
     * {@code BetaTest.testLogGammaMinusLogGammaSum()}, this private method is
     * accessed through reflection.
     *
     * @param double $a First argument.
     * @param double $b Second argument.
     * @return float the value of {@code log(Gamma(b) / Gamma(a + b))}.
     * @throws NumberIsTooSmallException if {@code a < 0.0} or {@code b < 10.0}.
     */
    private static function logGammaMinusLogGammaSum($a, $b)
    {
        if ($a < 0.0) {
            throw new NumberIsTooSmallException(a, 0.0, true);
        }
        if ($b < 10.0) {
            throw new NumberIsTooSmallException(b, 10.0, true);
        }

        /*
         * d = a + b - 0.5
         */
        $d = null;
        $w = null;
        if ($a <= $b) {
            $d = $b + ($a - 0.5);
            $w = $this->deltaMinusDeltaSum($a, $b);
        } else {
            $d = $a + ($b - 0.5);
            $w = $this->deltaMinusDeltaSum($b, $a);
        }

        $u = $d * log1p($a / $b);
        $v = $a * (log($b) - 1.0);

        return $u <= $v ? ($w - $u) - $v : ($w - $v) - $u;
    }

    /**
     * Returns the value of Δ(b) - Δ(a + b), with 0 ≤ a ≤ b and b ≥ 10. Based
     * on equations (26), (27) and (28) in Didonato and Morris (1992).
     *
     * @param $a First argument.
     * @param $b Second argument.
     * @return double The value of {@code Delta(b) - Delta(a + b)}
     * @throws OutOfRangeException if {@code a < 0} or {@code a > b}
     * @throws NumberIsTooSmallException if {@code b < 10}
     */
    private static function deltaMinusDeltaSum($, $b)
    {
        if (($a < 0) || ($a > $b)) {
            throw new OutOfRangeException($a, 0, $b);
        }
        if ($b < 10) {
            throw new NumberIsTooSmallException($b, 10, true);
        }

        $h = $a / $b;
        $p = $h / (1.0 + $h);
        $q = 1.0 / (1.0 + $h);
        $q2 = $q * $q;
        /*
         * s[i] = 1 + q + ... - q**(2 * i)
         */
        $s = array();
        $s[0] = 1.0;
        for ($i = 1; $i < count($s); $i++) {
            $s[$i] = 1.0 + ($q + $q2 * $s[$i - 1]);
        }
        /*
         * w = Delta(b) - Delta(a + b)
         */
        $sqrtT = 10.0 / $b;
        $t = $sqrtT * $sqrtT;
        $deltaLength = count($this->DELTA);
        $w = $this->DELTA[$deltaLength - 1] * $s[count($s)- 1];
        for ($i = $deltaLength - 2; $i >= 0; $i--) {
            $w = $t * $w + $this->DELTA[$i] * $s[$i];
        }
        return $w * $p / $b;
    }

    /**
     * Returns the value of Δ(p) + Δ(q) - Δ(p + q), with p, q ≥ 10. Based on
     * the <em>NSWC Library of Mathematics Subroutines</em> double precision
     * implementation, {@code DBCORR}. In
     * {@code BetaTest.testSumDeltaMinusDeltaSum()}, this private method is
     * accessed through reflection.
     *
     * @param $p First argument.
     * @param $q Second argument.
     * @return double The value of {@code Delta(p) + Delta(q) - Delta(p + q)}.
     * @throws NumberIsTooSmallException if {@code p < 10.0} or {@code q < 10.0}.
     */
    private static function sumDeltaMinusDeltaSum($p, $q)
    {
        if ($p < 10.0) {
            throw new NumberIsTooSmallException(p, 10.0, true);
        }
        if ($q < 10.0) {
            throw new NumberIsTooSmallException(q, 10.0, true);
        }

        $a = min($p, $q);
        $b = max($p, $q);
        $sqrtT = 10.0 / $a;
        $t = $sqrtT * $sqrtT;
        $deltaLength = count($this->DELTA);
        double z = $this->DELTA[$deltaLength - 1];
        for ($i = $deltaLength - 2; $i >= 0; $i--) {
            $z = $t * $z + $this->DELTA[$i];
        }
        return $z / $a + $this->deltaMinusDeltaSum($a, $b);
    }

    /**
     * Returns the value of log B(p, q) for 0 ≤ x ≤ 1 and p, q > 0. Based on the
     * <em>NSWC Library of Mathematics Subroutines</em> implementation,
     * {@code DBETLN}.
     *
     * @param float $p First argument.
     * @param float $q Second argument.
     * @return float The value of {@code log(Beta(p, q))}, {@code NaN} if
     * {@code p <= 0} or {@code q <= 0}.
     */
    public static function logBeta($p, $q)
    {
        if (NAN === $p || NAN === $q || ($p <= 0.0) || ($q <= 0.0)) {
            return NAN;
        }

        $a = min($p, $q);
        $b = max($p, $q);
        if ($a >= 10.0) {
            $w = $this->sumDeltaMinusDeltaSum($a, $b);
            $h = $a / $b;
            $c = $h / (1.0 + $h);
            $u = -($a - 0.5) * log($c);
            $v = $b * log1p($h);
            if ($u <= $v) {
                return (((-0.5 * FastMath.log(b) + HALF_LOG_TWO_PI) + w) - u) - v;
            }

            return (((-0.5 * FastMath.log(b) + HALF_LOG_TWO_PI) + w) - v) - u;
        }
        
        if ($a > 2.0) {
            if ($b > 1000.0) {
                $n = (int) floor($a - 1.0);
                $prod = 1.0;
                $ared = $a;
                for ($i = 0; $i < $n; $i++) {
                    $ared -= 1.0;
                    $prod *= $ared / (1.0 + $ared / $b);
                }
                return (log($prod) - $n * log($b)) +
                        (Gamma::logGamma($ared) +
                        $this->logGammaMinusLogGammaSum($ared, $b));
            }

            $prod1 = 1.0;
            $ared = $a;

            while ($ared > 2.0) {
                $ared -= 1.0;
                final double $h = $ared / $b;
                $prod1 *= $h / (1.0 + $h);
            }
            if ($b < 10.0) {
                $prod2 = 1.0;
                $bred = b;
                while ($bred > 2.0) {
                    $bred -= 1.0;
                    $prod2 *= $bred / ($ared + $bred);
                }
                return log($prod1) +
                       log($prod2) +
                       (Gamma::logGamma($ared) +
                       (Gamma::logGamma($bred) -
                       $this->logGammaSum($ared, $bred)));
            }

            return log($prod1) +
                   Gamma::logGamma($ared) +
                   $this->logGammaMinusLogGammaSum($ared, $b);
        } 
        
        if ($a >= 1.0) {
            if ($b > 2.0) {
                if ($b < 10.0) {
                    $prod = 1.0;
                    $bred = $b;
                    while ($bred > 2.0) {
                        $bred -= 1.0;
                        $prod *= $bred / ($a + $bred);
                    }
                    return log($prod) +
                           (Gamma::logGamma($a) +
                            (Gamma::logGamma($bred) -
                            $this->logGammaSum($a, $bred)));
                }
                return Gamma::logGamma($a) +
                    $this->logGammaMinusLogGammaSum($a, $b);
            }

            return Gamma::logGamma($a) +
                   Gamma::logGamma($b) -
                   $this->logGammaSum($a, $b);
        }

        if ($b >= 10.0) {
            return Gamma::logGamma($a) +
                   $this->logGammaMinusLogGammaSum($a, $b);
        }

        // The following command is the original NSWC implementation.
        // return Gamma.logGamma(a) +
        // (Gamma.logGamma(b) - Gamma.logGamma(a + b));
        // The following command turns out to be more accurate.
        return log(Gamma::gamma($a) * Gamma::gamma($b) /
                            Gamma::gamma($a + $b));
    }
}
