<?php

/**
 * Provides a generic means to evaluate continued fractions.  Subclasses simply
 * provided the a and b coefficients to evaluate the continued fraction.
 *
 * <p>
 * References:
 * <ul>
 * <li><a href="http://mathworld.wolfram.com/ContinuedFraction.html">
 * Continued Fraction</a></li>
 * </ul>
 * </p>
 *
 */
class ContinuedFraction 
{
    /** 
     * Maximum allowed numerical error.
     */
    const DEFAULT_EPSILON = 10e-9;

    /**
     * @var \Closure
     */
    private $aClosure;

    /**
     * @var \Closure
     */
    private $bClosure;

    public function __construct(\Closure $aClosure, \Closure $bClosure)
    {
        $this->aClosure = $aClosure;
        $this->bClosure = $bClosure;
    }

    /**
     * Evaluates the continued fraction at the value x.
     * <p>
     * The implementation of this method is based on the modified Lentz algorithm as described
     * on page 18 ff. in:
     * <ul>
     *   <li>
     *   I. J. Thompson,  A. R. Barnett. "Coulomb and Bessel Functions of Complex Arguments and Order."
     *   <a target="_blank" href="http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf">
     *   http://www.fresco.org.uk/papers/Thompson-JCP64p490.pdf</a>
     *   </li>
     * </ul>
     * <b>Note:</b> the implementation uses the terms a<sub>i</sub> and b<sub>i</sub> as defined in
     * <a href="http://mathworld.wolfram.com/ContinuedFraction.html">Continued Fraction @ MathWorld</a>.
     * </p>
     *
     * @param float $x the evaluation point.
     * @param float $epsilon maximum error allowed.
     * @param float $maxIterations maximum number of convergents
     * @return float The value of the continued fraction evaluated at x.
     * @throws ConvergenceException if the algorithm fails to converge.
     * @throws MaxCountExceededException if maximal number of iterations is reached
     */
    public function evaluate($x, $epsilon = self::DEFAULT_EPSILON, $maxIterations = PHP_INT_MAX)
    {
        // should be 1e-50
        //$small = 0.0000000000000000000000000000000000000000000000001;

        $hPrev = $this->aClosure(0, $x);

        // use the value of small as epsilon criteria for zero checks
        //if (Precision.equals(hPrev, 0.0, small)) {
            //hPrev = small;
        //}

        $n = 1;
        $dPrev = 0.0;
        $cPrev = $hPrev;
        $hN = $hPrev;

        while ($n < $maxIterations) {
            $a = $this->aClosure($n, $x);
            $b = $this->bClosure($n, $x);

            $dN = $a + $b * $dPrev;
            //if (Precision.equals(dN, 0.0, small)) {
                //dN = small;
            //}
            $cN = $a + $b / $cPrev;
            //if (Precision.equals(cN, 0.0, small)) {
                //cN = small;
            //}

            $dN = 1 / $dN;
            $deltaN = $cN * $dN;
            $hN = $hPrev * $deltaN;

            // really not sure how this applies
            if (INF === $hN) {
                throw new ConvergenceException('Encountered infinity divergence..');
            }

            // really not sure how this applies
            if (NAN === $hN) {
                throw new ConvergenceException('Encountered non-numeric result..');
            }

            if (abs($deltaN - 1.0) < $epsilon) {
                break;
            }

            $dPrev = $dN;
            $cPrev = $cN;
            $hPrev = $hN;
            $n++;
        }

        if ($n >= $maxIterations) {
            throw new MaxCountExceededException(sprintf(
                'Maximum number of iterations (%s) exceeded for evaluation point %d',
                $maxIterations, $x
            ));
        }

        return $hN;
    }
}
