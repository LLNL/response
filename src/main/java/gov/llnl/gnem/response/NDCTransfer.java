/*-
 * #%L
 * Seismic Response Processing Module
 *  LLNL-CODE-856351
 *  This work was performed under the auspices of the U.S. Department of Energy
 *  by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.
 * %%
 * Copyright (C) 2023 Lawrence Livermore National Laboratory
 * %%
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * #L%
 */

package gov.llnl.gnem.response;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.math3.complex.Complex;

import com.isti.jevalresp.ResponseUnits;

/**
 *
 * @author dodge1
 */
public class NDCTransfer {

    private static final double TWOPI = 2 * Math.PI;
 
    public TransferData getFromTransferFunction(int nsamp, double samprate, double time, ResponseMetaData metadata) throws IOException {
        int nfft = TransferFunctionUtils.next2(nsamp);
        int nfreq = nfft / 2 + 1;

        double[] xre = new double[nfreq];
        double[] xim = new double[nfreq];

        double delta = 1.0 / samprate;
        double delfrq = 1.0 / (nfft * delta);
        ResponseType fromType = metadata.getRsptype();

        computeTransferFunctionFromFile(nfreq, delfrq, time, xre, xim, fromType, metadata.getFilename());
        NDCUnitsParser parser = new NDCUnitsParser();
        ResponseUnits ru = parser.getResponseUnits(new File(metadata.getFilename()));
        return TransferFunctionUtils.packageTransferFunctionSamples(delfrq, xre, xim, ru, metadata);
    }

    public void computeTransferFunctionFromFile(int nfreq, double delfrq, double epoch, double[] xre, double[] xim, ResponseType type, String filename) throws IOException {
        /*
         * Branching routine to apply the individual instrument responses.
         */
        for (int idx = 0; idx < nfreq; idx++) {
            xre[idx] = 1.0e0;
            xim[idx] = 0.0e0;
        }
        switch (type) {
            case CSS:
            case PAZ:
            case FAP:
            case PAZFIR:
                transfer(filename, delfrq, nfreq, xre, xim);
                break;
            default:
                throw new IllegalStateException("Unhandled type: " + type);
        }

    }

    public static void transfer(String filename, double dt, int nfr, double[] xre, double[] xim) throws FileNotFoundException {
        double startFrequency = 0.0;
        double endFrequency = dt * (nfr - 1);

        Polar[] cascade = new Polar[nfr];
        for (int i = 0; i < cascade.length; i++) {
            cascade[i] = new Polar(0.0, 1.0);
        }

        Scanner sc = new Scanner(new File(filename));
        Double isr = null;
        while (sc.hasNextLine()) {
            String line = sc.nextLine();
            if (line == null || line.trim().isEmpty()) {
                continue;
            }
            if (line.charAt(0) != '#') {
                if (line.contains("theoretical") || line.contains("measured")) {
                    Polar[] result = new Polar[0];
                    if (line.contains("paz")) {
                        result = processNDCpaz(sc, result, nfr, startFrequency, endFrequency);
                    } else if (line.contains("PAZ2")) {
                        result = processIDCpaz(line, sc, result, nfr, startFrequency, endFrequency);
                    } else if (line.contains("fap")) {
                        Fap[] faps = readFap(sc);
                        result = doFap(nfr, startFrequency, endFrequency, faps);
                    } else if (line.contains("FAP2")) {
                        Fap[] faps = readIDCFap(line, sc);
                        result = doFap(nfr, startFrequency, endFrequency, faps);
                    } else if (line.contains("fir")) {
                        Fir fir = readFir(sc);
                        result = doFir(nfr, startFrequency, endFrequency, fir);
                    } else if (line.contains("DIG2")) {
                        isr = processDig2Line(line);
                        continue;
                    } else if (line.contains("FIR2")) {
                        if (isr == null) {
                            throw new IllegalStateException("Encountered FIR2 line not preceeded by DIG2 line!");
                        }
                        Fir fir = readIDCFir(line, isr, sc);
                        result = doFir(nfr, startFrequency, endFrequency, fir);

                    }

                    /*
                     * Cascade individual group responses
                     */
                    for (int j = 0; j < result.length; j++) {
                        cascade[j].a *= result[j].a;
                        cascade[j].p += result[j].p;
                    }
                }
            }
        }

        for (int j = 0; j < nfr; j++) {
            Complex c = tocmplx(cascade[j]);
            xre[j] = c.getReal();
            xim[j] = c.getImaginary();
        }
    }

    private static double readNormFactor(Scanner sc) {
        double normFactor = sc.nextDouble();
        sc.nextLine();
        return normFactor;
    }

    private static DComplex[] readDComplex(Scanner sc) {
        String line = sc.nextLine();

        int linesLeft = getFirstInt(line);
        DComplex[] results = new DComplex[linesLeft];
        while (linesLeft > 0) {
            StringTokenizer tokenizer = new StringTokenizer(sc.nextLine());

            double zReal = Double.parseDouble(tokenizer.nextToken());
            double zImag = Double.parseDouble(tokenizer.nextToken());
            double eReal = Double.parseDouble(tokenizer.nextToken());
            double eImag = Double.parseDouble(tokenizer.nextToken());

            DComplex dc = new DComplex(new Complex(zReal, zImag), new Complex(eReal, eImag));
            results[results.length - linesLeft] = dc;

            linesLeft--;
        }

        return results;
    }

    private static DComplex[] readDComplexIDC(Scanner sc, int nValues) {
        DComplex[] results = new DComplex[nValues];
        for (int j = 0; j < nValues; ++j) {
            StringTokenizer tokenizer = new StringTokenizer(sc.nextLine());

            double zReal = Double.parseDouble(tokenizer.nextToken());
            double zImag = Double.parseDouble(tokenizer.nextToken());
            double eReal = 0.0;
            double eImag = 0.0;
            DComplex dc = new DComplex(new Complex(zReal, zImag), new Complex(eReal, eImag));
            results[j] = dc;
        }

        return results;
    }

    private static Fap[] readFap(Scanner sc) {
        String line = sc.nextLine();

        int linesLeft = getFirstInt(line);
        Fap[] results = new Fap[linesLeft];
        while (linesLeft > 0) {
            StringTokenizer tokenizer = new StringTokenizer(sc.nextLine());

            double f = Double.parseDouble(tokenizer.nextToken());
            double a = Double.parseDouble(tokenizer.nextToken());
            double p = Double.parseDouble(tokenizer.nextToken());
            double ae = Double.parseDouble(tokenizer.nextToken());
            double pe = Double.parseDouble(tokenizer.nextToken());

            p *= TWOPI / 360.0;

            Fap fap = new Fap(f, a, p, ae, pe);
            results[results.length - linesLeft] = fap;

            linesLeft--;
        }

        return results;
    }

    private static Fir readFir(Scanner sc) {
        double isr = sc.nextDouble();
        Fir fir = new Fir(isr);
        sc.nextLine();

        String line = sc.nextLine();

        int linesLeft = getFirstInt(line);
        int numCoefficients = linesLeft;
        fir.setNumerator(linesLeft);
        double sum = 0;
        while (linesLeft > 0) {
            StringTokenizer tokenizer = new StringTokenizer(sc.nextLine());

            double nu = Double.parseDouble(tokenizer.nextToken());
            double nue = Double.parseDouble(tokenizer.nextToken());

            int index = fir.nnc - linesLeft;
            fir.nu[index] = nu;
            fir.nue[index] = nue;

            linesLeft--;
            sum += nu;
        }

        line = sc.nextLine();
        linesLeft = Integer.parseInt(line.trim());
        fir.setDenominator(linesLeft);
        while (linesLeft > 0) {
            StringTokenizer tokenizer = new StringTokenizer(sc.nextLine());

            double de = Double.parseDouble(tokenizer.nextToken());
            double dee = Double.parseDouble(tokenizer.nextToken());

            int index = fir.ndc - linesLeft;
            fir.de[index] = de;
            fir.dee[index] = dee;

            linesLeft--;
        }
        String msg = String.format("Created FIR with ISR = %f and %d coeficients. Sum of coefficients is %f",
                isr, numCoefficients, sum);
        Logger.getLogger("NDCTransfer").log(Level.FINE, msg);

        return fir;
    }

    private static Fir readIDCFir(String line, Double isr, Scanner sc) {

        String[] tokens = line.split("\\s+");
        if (tokens.length < 7) {
            throw new IllegalStateException("Expected FIR2 line to have at least 7 tokens, but line is: " + line + "!");
        }
        Fir fir = new Fir(isr);
        int numCoefficients = Integer.parseInt(tokens[6]);
        fir.setNumerator(numCoefficients);
        int index = 0;
        double sum = 0;
        while (index < numCoefficients) {
            line = sc.nextLine();
            tokens = line.split("\\s+");
            for (String token : tokens) {
                if (!token.isEmpty()) {
                    double nu = Double.parseDouble(token);
                    double nue = 0.0;
                    fir.nu[index] = nu;
                    fir.nue[index++] = nue;
                    sum += nu;
                }
            }
        }

        fir.setDenominator(0);
        String msg = String.format("Created FIR with ISR = %f and %d coeficients. Sum of coefficients is %f",
                isr, numCoefficients, sum);
        Logger.getLogger("NDCTransfer").log(Level.FINE, msg);
        return fir;
    }

    /*
     * Compute response using poles and zeros
     */
    private static Polar[] doPaz(int nfr, double start_fr, double end_fr, DComplex[] poles, DComplex[] zeros) {
        double amp, phase, delta_f;

        if (nfr == 1) /*
                       * avoid 0 in denominator
         */ {
            delta_f = 1.0;
            /*
                            * value irrelevant, loop ends first
             */
        } else {
            /*
             * constant linear spacing
             */
            delta_f = (end_fr - start_fr) / (nfr - 1.0);
        }

        Polar[] result = new Polar[nfr];
        double delta = 1.0;
        for (int j = 0; j < nfr; j++, delta *= delta_f) {
            double omega = TWOPI * (start_fr + j * delta_f);
            result[j] = new Polar(0.0, 1.0);

            for (int i = 0; i < zeros.length; i++) {
                Polar polar = topolar(-zeros[i].z.getReal(), (omega - zeros[i].z.getImaginary()));
                result[j].a *= polar.a;
                result[j].p += polar.p;
            }

            for (int i = 0; i < poles.length; i++) {
                Polar polar = topolar(-poles[i].z.getReal(), (omega - poles[i].z.getImaginary()));
                if (polar.a != 0.0) {
                    result[j].a /= polar.a;
                }
                result[j].p -= polar.p;
            }
        }

        return result;
    }

    private static Polar[] processNDCpaz(Scanner sc, Polar[] result, int nfr, double start_fr, double end_fr) {
        double normFactor = readNormFactor(sc);
        DComplex[] poles = readDComplex(sc);
        DComplex[] zeros = readDComplex(sc);
        result = doPaz(nfr, start_fr, end_fr, poles, zeros);
        for (int j = 0; j < nfr; j++) {
            result[j].a *= normFactor;
        }
        return result;
    }

    private static Polar[] processIDCpaz(String line, Scanner sc, Polar[] result, int nfr, double startFrequency, double endFrequency) {
        String[] tokens = line.split("\\s+");
        if (tokens.length < 8) {
            throw new IllegalStateException("Expected PAZ2 line to have at least 8 tokens, but line is: " + line + "!");
        }
        double normFactor = Double.parseDouble(tokens[3]);
        int numPoles = Integer.parseInt(tokens[6]);
        int numZeros = Integer.parseInt(tokens[7]);
        DComplex[] poles = readDComplexIDC(sc, numPoles);
        DComplex[] zeros = readDComplexIDC(sc, numZeros);
        result = doPaz(nfr, startFrequency, endFrequency, poles, zeros);
        for (int j = 0; j < nfr; j++) {
            result[j].a *= normFactor;
        }
        return result;
    }

    private static Double processDig2Line(String line) {
        String[] tokens = line.split("\\s+");
        if (tokens.length < 4) {
            throw new IllegalStateException("Expected DIG2 line to have at least 4 tokens, but line is: " + line + "!");
        }
        return Double.valueOf(tokens[3]);
    }

    private static Fap[] readIDCFap(String line, Scanner sc) {
        String[] tokens = line.split("\\s+");
        if (tokens.length < 6) {
            throw new IllegalStateException("Expected FAP2 line to have at least 6 tokens, but line is: " + line + "!");
        }
        int nfrequencies = Integer.parseInt(tokens[5]);

        Fap[] results = new Fap[nfrequencies];
        int frequenciesRead = 0;
        while (frequenciesRead < nfrequencies) {
            line = sc.nextLine();
            tokens = line.trim().split("\\s+");
            if (tokens.length == 3) {
                double f = Double.parseDouble(tokens[0]);
                double a = Double.parseDouble(tokens[1]);
                double p = Double.parseDouble(tokens[2]);
                double ae = 0.0;
                double pe = 0.0;
                p *= TWOPI / 360.0;
                Fap fap = new Fap(f, a, p, ae, pe);
                results[frequenciesRead++] = fap;
            }
        }

        return results;

    }

    private static Polar[] doFap(int numRequestedFrequencies, double requestedStartFrequency, double requestedEndFrequency, Fap[] faps) {
        int order;
        double freq;
        /*
                      * frequency or log(f) if log_flag
         */
        double ftmp1;
        double ftmp2;
        double requestedDeltaF;

        if (numRequestedFrequencies == 1) {
            /*
             * avoid 0 in denominator
             */
            requestedDeltaF = 1.0;
        } else {
            requestedDeltaF = (requestedEndFrequency - requestedStartFrequency) / (numRequestedFrequencies - 1.0);
        }

        int numFaps = faps.length;
        double maxFreqInFapFile = faps[numFaps - 1].f;
        if (requestedEndFrequency > maxFreqInFapFile) {
            Fap lastInFile = faps[numFaps - 1];
            Fap[] extended = new Fap[numFaps + 1];
            for (int j = 0; j < numFaps; ++j) {
                extended[j] = faps[j];
            }
            extended[numFaps] = new Fap(requestedEndFrequency, lastInFile.a, lastInFile.p, lastInFile.ae, lastInFile.pe);
            faps = extended;
        }
        /*
         * Lagrange interpolation works better for log function
         */
        double[] f = new double[faps.length];
        double[] a = new double[faps.length];
        for (int i = 0; i < faps.length; i++) {
            if (faps[i].f <= 1.0e-20) {
                faps[i].f = 1.0e-20; //To avoid log(0)
            }
            if (faps[i].a <= 1.0e-20) {
                faps[i].a = 1.0e-20;//To avoid log(0)
            }

            f[i] = Math.log10(faps[i].f);
            a[i] = Math.log10(faps[i].a);
        }

        Polar[] result = new Polar[numRequestedFrequencies];
        for (int j = 0; j < numRequestedFrequencies; j++) {
            ftmp2 = requestedStartFrequency + j * requestedDeltaF;
            if (ftmp2 <= 1.0e-20) {
                ftmp2 = 1.0e-20;//To avoid log(0)
            }

            freq = Math.log10(ftmp2);
            result[j] = new Polar(0.0, 1.0);
            int i;
            for (i = 0; i < faps.length; i++) {
                if (freq < f[i]) {
                    break;
                }
            }

            i -= 2;
            order = 4;
            if (i < 0) {
                i = 0;
                order = 2;
            } else if (i > faps.length - 4) {
                i = faps.length - 2;
                order = 2;
            }

            double[] p = new double[faps.length];
            for (int k = 0; k < faps.length; k++) {
                p[k] = faps[k].p;
            }

            ftmp1 = lagrange(a, f, i, order, freq);
            result[j].a = ftmp1;
            ftmp1 = lagrange(p, f, i, order, freq);
            result[j].p = ftmp1;
            result[j].a = Math.pow(10.0, result[j].a);
        }
        return result;
    }

    private static Polar[] doFir(int nfr, double start_fr, double end_fr, Fir firs) {
        int flag = 1;
        /*
                       * 1= forward transform, -1 reverse
         */

 /*
         * Initialize for error handling, one return at bailout
         */

 /*
         * frequency amplitude phase
         */
        int n = 513;
        Fap[] faps = new Fap[n];
        double[] xr = new double[n * 2];

        for (int i = 0; i < n; i++) {
            faps[i] = new Fap(0.0, 0.0, 0.0, 0.0, 0.0);
        }

        /*
         * Set up data in a large array
         */
        for (int i = 0; i < firs.nnc; i++) {
            xr[i] = firs.nu[i];
        }
        for (int i = firs.nnc; i < faps.length * 2; i++) {
            xr[i] = 0.0;
        }

        /*
         * Set up and call fft for numerator coefficients
         */
        double df = 1.0;

        /*
         * Find greatest power of 2 that is less than twice faps.n
         */
        int nexp;
        for (nexp = 1; Math.pow(2.0, nexp) < 2 * faps.length; nexp++) {
        }
        nexp--;

        /*
         * note df no longer passed since no scaling is to be done, ganz 6/11
         */
        odfftr(nexp, xr, flag);

        /*
         * Compute frequency, amplitude, and phase for numerator (fir)
         */
        df = firs.isr / ((faps.length - 1) * 2);
        for (int i = 0; i < faps.length; i++) {
            faps[i].a = Math.hypot(xr[i * 2], xr[i * 2 + 1]);
            if (xr[i * 2] == 0 && xr[i * 2 + 1] == 0) {
                faps[i].p = 0;
            } else {
                faps[i].p = Math.atan2(xr[i * 2], xr[i * 2 + 1]);
            }
            faps[i].f = i * df;
        }

        if (firs.ndc > 0) {
            /*
             * Set up data in a large array
             */
            for (int i = 0; i < firs.ndc; i++) {
                xr[i] = firs.de[i];
            }
            for (int i = firs.ndc; i < faps.length * 2; i++) {
                xr[i] = 0.0;
            }

            /*
             * Set up and call fft for numerator coefficients
             */
            df = 1.0;
            /*
             * df for scaling, no longer passed since no scaling was ever done
             */
            odfftr(nexp, xr, flag);

            /*
             * Compute frequency, amplitude, and phase for denom. (iir)
             */
            for (int i = 0; i < faps.length; i++) {
                faps[i].a /= Math.hypot((1.0 - xr[i * 2]), xr[i * 2 + 1]);
                if ((1.0 - xr[i * 2]) == 0 && xr[i * 2 + 1] == 0) {
                    faps[i].p -= 0;
                } else {
                    faps[i].p -= Math.atan2((1.0 - xr[i * 2]), xr[i * 2 + 1]);
                }
            }
        }

        /*
         * Make the phase a smooth function, i.e. get rid of wraps \/\
         */
 /*
         * Find the index with the derivative closest to zero
         */
        // TODO this is a bug in the C version of the code, change the type of min
        // back to double after testing is complete
        // double min = TWOPI;
        int min = (int) TWOPI;
        int minindex = 0;
        double[] deriv = new double[faps.length];
        for (int i = 0; i < faps.length - 1; i++) {
            deriv[i] = faps[i + 1].p - faps[i].p;
            if (Math.abs(deriv[i]) < min) {
                min = (int) Math.abs(deriv[i]);
                minindex = i;
            }
        }

        /*
         * Correct phases with high indicies
         */
        for (int i = minindex; i < faps.length - 2; i++) {
            if (sign(deriv[i]) != sign(deriv[i + 1])) {
                for (int j = i + 2; j < faps.length; j++) {
                    faps[j].p = faps[j].p + sign(deriv[i]) * TWOPI;
                }
                deriv[i + 1] = faps[i + 2].p - faps[i + 1].p;
            }
        }

        /*
         * Correct phases with low indicies
         */
        for (int i = minindex; i > 0; i--) {
            if (sign(deriv[i]) != sign(deriv[i - 1])) {
                for (int j = i - 1; j >= 0; j--) {
                    faps[j].p = faps[j].p + sign(deriv[i]) * TWOPI;
                }
                deriv[i - 1] = faps[i].p - faps[i - 1].p;
            }
        }

        /*
         * Call fap to interpolate for proper values
         */
        return doFap(nfr, start_fr, end_fr, faps);
    }

    private static Polar topolar(double real, double imag) {
        double amp = Math.sqrt(real * real + imag * imag);
        double phase;
        if (imag == 0.0 && real == 0.0) {
            phase = 0.0;
        } else {
            phase = Math.atan2(imag, real);
        }

        return new Polar(phase, amp);
    }

    private static double lagrange(double[] f, double[] xi, int offset, int n, double x) {
        /*
         * evenly spaced , in order
         */
        int i, k;
        double prod;

        double fx = 0.0;
        for (k = 0; k < n; k++) {
            prod = 1.0;
            for (i = 0; i < n; i++) {
                if (i != k) {
                    prod = prod * (x - xi[i + offset]) / (xi[k + offset] - xi[i + offset]);
                }
            }
            fx = fx + prod * f[k + offset];
        }
        return fx;
    }

    private static void odfftr(int nexp, double[] xr, int flag) {
        /*
         * complex c;
         */
        double c;
        int i, j, k, m, p, n1;
        int Ls, ks, ms, jm, dk;
        double wr, wi, tr, ti;

        double[] xi;
        int n2;
        int npts, n, jj;

        /*
         * need to build imaginary array
         */
        npts = 1;
        n2 = 1;

        npts = 1;

        i = nexp;
        while (i-- > 0) {
            npts <<= 1;
        }

        n = npts;
        xi = new double[npts + 2];

        if (flag < 0) {
            for (i = 0; i <= npts / 2; i++) {
                k = i * 2;
                xi[i] = xr[k + 1];
                xr[i] = xr[k];
            }
            for (i = 0; i <= npts / 2; i++) {
                k = i * 2;
                xi[npts - i] = -xi[i];
                xr[npts - i] = xr[i];
            }

        } else {
            for (i = 0; i < npts; i++) {
                xi[i] = 0.0f;
            }

        }

        n1 = n / n2;
        /*
         * do bit reversal permutation
         */
        for (k = 0; k < n1; ++k) {
            /*
                                    * This is algorithms 1.5.1 and 1.5.2.
             */
            j = 0;
            m = k;
            p = 1;
            /*
             * p = 2^q, q used in the book
             */
            while (p < n1) {
                j = 2 * j + (m & 1);
                m >>= 1;
                p <<= 1;
            }

            assert (p == n1);
            /*
             * make sure n1 is a power of two
             */
            if (j > k) {
                for (i = 0; i < n2; ++i) {
                    /*
                     * swap k <-> j row
                     */
                    c = xr[k * n2 + i];
                    /*
                     * for all columns
                     */
                    xr[k * n2 + i] = xr[j * n2 + i];
                    xr[j * n2 + i] = c;

                    c = xi[k * n2 + i];
                    /*
                     * for all columns
                     */
                    xi[k * n2 + i] = xi[j * n2 + i];
                    xi[j * n2 + i] = c;

                }
            }
        }

        /*
         * This is (3.1.7), page 124
         */
        p = 1;
        while (p < n1) {
            Ls = p;
            p <<= 1;
            jm = 0;
            /*
             * jm is j*n2
             */
            dk = p * n2;
            for (j = 0; j < Ls; ++j) {
                wr = Math.cos(Math.PI * j / Ls);
                /*
                 * real and imaginary part
                 */
                wi = -flag * Math.sin(Math.PI * j / Ls);
                /*
                 * of the omega
                 */
                for (k = jm; k < n; k += dk) {
                    /*
                     * "butterfly"
                     */
                    ks = k + Ls * n2;
                    for (i = 0; i < n2; ++i) {
                        /*
                         * for each row
                         */
                        m = k + i;
                        ms = ks + i;
                        tr = wr * xr[ms] - wi * xi[ms];
                        ti = wr * xi[ms] + wi * xr[ms];
                        xr[ms] = xr[m] - tr;
                        xi[ms] = xi[m] - ti;
                        xr[m] += tr;
                        xi[m] += ti;
                    }
                }
                jm += n2;
            }
        }

        /*
         * now combine the real/imaginary parts back into xr
         */
        if (flag > 0) {
            for (i = npts / 2; i >= 0; i--) {
                k = i * 2;
                xr[k] = xr[i];
                xr[k + 1] = xi[i];
            }
        }
    }

    private static int sign(double a) {
        return a < 0.0 ? -1 : 1;
    }

    private static Complex tocmplx(Polar polar) {
        return new Complex(polar.a * Math.cos(polar.p), polar.a * Math.sin(polar.p));
    }

    private static int getFirstInt(String line) {
        Pattern intPattern = Pattern.compile("\\d+");
        Matcher matcher = intPattern.matcher(line);
        matcher.find();
        return Integer.parseInt(matcher.group());
    }

    private static class DComplex {

        public final Complex z;
        public final Complex e;

        public DComplex(Complex z, Complex e) {
            this.z = z;
            this.e = e;
        }
    }

    private static class Fap {

        public double f;
        public double a;
        public double p;
        public double ae;
        public double pe;

        public Fap(double f, double a, double p, double ae, double pe) {
            this.f = f;
            this.a = a;
            this.p = p;
            this.ae = ae;
            this.pe = pe;
        }
    }

    private static class Fir {

        private double isr;
        public int nnc;
        public int ndc;
        public double[] nu;
        public double[] nue;
        public double[] de;
        public double[] dee;

        public Fir(double isr) {
            this.isr = isr;
        }

        public void setNumerator(int nnc) {
            this.nnc = nnc;
            nu = new double[nnc];
            nue = new double[nnc];
        }

        public void setDenominator(int ndc) {
            this.ndc = ndc;
            de = new double[ndc];
            dee = new double[ndc];
        }

    }

    private static class Polar {

        public double p;
        public double a;

        public Polar(double p, double a) {
            this.p = p;
            this.a = a;
        }

        @Override
        public String toString() {
            return "Polar{" + "p=" + p + ", a=" + a + '}';
        }

    }

}
