//
//  ATSMath.h
//  AllTrailsShared
//
//  Created by Theodore Calmes on 5/20/13.
//  Copyright (c) 2013 AllTrails. All rights reserved.
//

#import <Accelerate/Accelerate.h>

static double ATSDegreesToRadian(double degrees)
{
    return (degrees / 180.0) * M_PI;
}

static double ATSRadiansToDegrees(double radians)
{
    return (radians / M_PI) * 180.0;
}

static double ATSHaversineDistance(double lat1, double lat2, double lon1, double lon2)
{
    double R = 6378100;
    double dLat = ATSDegreesToRadian(lat2 - lat1);
    double dLon = ATSDegreesToRadian(lon2 - lon1);
    double l1 = ATSDegreesToRadian(lat1);
    double l2 = ATSDegreesToRadian(lat2);

    double a = pow(sin(dLat/2.0), 2.0) + pow(sin(dLon/2.0), 2.0) * cos(l1) * cos(l2);
    double c = 2 * atan(sqrt(a) / sqrt(1-a));

    return R * c;
}

static CGFloat ATSDistanceBetweenPoints(CGPoint a, CGPoint b)
{
    return sqrtf(powf(a.x - b.x, 2) + powf(a.y - b.y, 2));
}

static CGFloat ATSAngleBetweenPoints(CGPoint a, CGPoint b)
{
    return atan2f(b.y - a.y, b.x - a.x);
}

static void printMatrix(float *m, int rows, int cols)
{
    printf("\nPrinting Matrix\n");
    for (NSInteger i = 0; i < rows; i++) {
        for (NSInteger j = 0; j < cols; j++) {
            printf(" %e ", m[i * cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

static float *ATSMatrixMakeWithDiagonals(float *diagonals, int count)
{
    float *matrix = (float *)calloc(count * count, sizeof(float));
    for (NSInteger i = 0; i < count; i++) {
        for (NSInteger j = 0; j < count; j++) {
            if (i == j) {
                matrix[i * count + j] = diagonals[i];
            }
        }
    }

    return matrix;
}

static double *ATSMatrixMakeWithDiagonal(double diagonal, int count)
{
    double *matrix = (double *)calloc(count * count, sizeof(double));
    for (NSInteger i = 0; i < count; i++) {
        matrix[i * count + i] = diagonal;
    }

    return matrix;
}

static float ATSTotal(float *values, int count)
{
    float sum = 0;
    vDSP_sve(values, 1, &sum, count);

    return sum;
}

static float ATSMinimum(float *values, int count)
{
    float min = FLT_MAX;
    for (int i = 0; i < count; i++) {
        if (min > values[i]) {
            min = values[i];
        }
    }
    return min;
}

static float ATSMaximum(float *values, int count)
{
    float max = -FLT_MAX;
    for (int i = 0; i < count; i++) {
        if (max < values[i]) {
            max = values[i];
        }
    }
    return max;
}

static float ATSAverage(float *values, int count)
{
    float avg = 0;
    vDSP_meanv(values, 1, &avg, count);

    return avg;
}

static float ATSStandardDeviation(float *values, int count)
{
    if (count <= 1) {
        return 0;
    }
    
    float avg = -ATSAverage(values, count);
    float *temp = malloc(count * sizeof(float));

    float deviation = 0;

    vDSP_vsadd(values, 1, &avg, temp, 1, count);
    vDSP_svesq(temp, 1, &deviation, count);

    deviation = powf(deviation / (count - 1), 0.5);

    free(temp);

    return deviation;
}

static float ATSVariance(float *values, int count)
{
    return powf(ATSStandardDeviation(values, count), 2);
}

static float *ATSMatrixInverse(float *matrix, int count)
{
    long N = (long)count;
    float *A = malloc(count * count * sizeof(float));
    for (int i = 0; i < count * count; i++) {
        A[i] = matrix[i];
    }

    long ERR;
    long LDWORK = -1;
    float *WORK = calloc(1, sizeof(float));
    long NRHS = count;
    long LDA = count;
    long LDB = count;

    float *sing = calloc(count * count, sizeof(float));
    float RCOND = -1.0f;
    long IRANK = -1;

    // Identity Matrix
    float *B = calloc(count * count, sizeof(float));
    for(int i=0; i<count; i++)  {
        B[i*count + i] = 1.0;
    }

    // sgelss solves for the least squares equation using singular value decomposition. By passing b = Identity(dimmentions) to |Ax - b| we get back the inverse of A.
    sgelss_(&N, &N, &NRHS, A, &LDA, B, &LDB, sing, &RCOND, &IRANK, WORK, &LDWORK, &ERR);

    LDWORK = WORK[0];
    free(WORK);
    WORK = calloc(LDWORK, sizeof(float));

    sgelss_(&N, &N, &NRHS, A, &LDA, B, &LDB, sing, &RCOND, &IRANK, WORK, &LDWORK, &ERR);

    free(A);
    free(sing);
    free(WORK);

    return B;
}

static double *ATSMatrixInverseD(double *matrix, int count)
{
    long N = (long)count;
    double *A = malloc(count * count * sizeof(double));
    for (int i = 0; i < count * count; i++) {
        A[i] = matrix[i];
    }

    long ERR;
    long LDWORK = -1;
    double *WORK = calloc(1, sizeof(double));
    long NRHS = count;
    long LDA = count;
    long LDB = count;

    double *sing = calloc(count * count, sizeof(double));
    double RCOND = -1.0f;
    long IRANK = -1;

    // Identity Matrix
    double *B = calloc(count * count, sizeof(double));
    for(int i=0; i<count; i++)  {
        B[i*count + i] = 1.0;
    }

    // sgelss solves for the least squares equation using singular value decomposition. By passing b = Identity(dimmentions) to |Ax - b| we get back the inverse of A.
    dgelss_(&N, &N, &NRHS, A, &LDA, B, &LDB, sing, &RCOND, &IRANK, WORK, &LDWORK, &ERR);

    LDWORK = WORK[0];
    free(WORK);
    WORK = calloc(LDWORK, sizeof(double));

    dgelss_(&N, &N, &NRHS, A, &LDA, B, &LDB, sing, &RCOND, &IRANK, WORK, &LDWORK, &ERR);

    free(A);
    free(sing);
    free(WORK);

    return B;
}

static float *ATSMeansVector(float *values, int dimensions, int count)
{
    float *transposed = malloc(dimensions * count * sizeof(float));
    vDSP_mtrans(values, 1, transposed, 1, dimensions, count);

    float *ones = (float *)calloc(dimensions, sizeof(float));
    float one = 1.0;
    vDSP_vsadd(ones, 1, &one, ones, 1, dimensions);

    float *means = malloc(dimensions * sizeof(float));
    float N = (float)count;
    vDSP_mmul(transposed, 1, ones, 1, means, 1, dimensions, 1, count);
    vDSP_vsdiv(means, 1, &N, means, 1, dimensions);

    free(transposed);
    free(ones);

    return means;
}

static float *ATSCovarianceMatrix(float *values, int dimensions, int count)
{
    float *transposed = malloc(dimensions * count * sizeof(float));
    vDSP_mtrans(values, 1, transposed, 1, dimensions, count);

    float *ones = (float *)calloc(dimensions, sizeof(float));
    float one = 1.0;
    vDSP_vsadd(ones, 1, &one, ones, 1, dimensions);

    float *means = malloc(dimensions * sizeof(float));
    float N = (float)count;
    vDSP_mmul(transposed, 1, ones, 1, means, 1, dimensions, 1, count);
    vDSP_vsdiv(means, 1, &N, means, 1, dimensions);

    float *meansMatrix = malloc(dimensions * count * sizeof(float));
    vDSP_mmul(ones, 1, means, 1, meansMatrix, 1, count, dimensions, 1);

    float *meansDifference = malloc(dimensions * count * sizeof(float));
    vDSP_vsub(values, 1, meansMatrix, 1, meansDifference, 1, dimensions * count);

    float *meansDifferenceTransposed = malloc(dimensions * count * sizeof(float));
    vDSP_mtrans(meansDifference, 1, meansDifferenceTransposed, 1, dimensions, count);

    float *covarianceMatrix = malloc(dimensions * dimensions * sizeof(float));
    vDSP_mmul(meansDifferenceTransposed, 1, meansDifference, 1, covarianceMatrix, 1,dimensions, dimensions, count);

    float adjustedPolulation = (count - 1);

    vDSP_vsdiv(covarianceMatrix, 1, &adjustedPolulation, covarianceMatrix, 1, dimensions * dimensions);

    free(transposed);
    free(ones);
    free(means);
    free(meansMatrix);
    free(meansDifference);
    free(meansDifferenceTransposed);

    return covarianceMatrix;
}

static float ATSMahalanobisDistance(float *x, float *values, int dimensions, int count)
{
    float *means = ATSMeansVector(values, dimensions, count);
    float *covarianceMatrix = ATSCovarianceMatrix(values, dimensions, count);

    float *meanDifference = malloc(dimensions * sizeof(float));
    vDSP_vsub(x, 1, means, 1, meanDifference, 1, dimensions);

    float *covarianceMatrixInverse = ATSMatrixInverse(covarianceMatrix, dimensions);

    float *workspace = malloc(dimensions * sizeof(float));
    vDSP_mmul(covarianceMatrixInverse, 1, meanDifference, 1, workspace, 1, dimensions, 1, dimensions);

    float distance = 0;
    vDSP_dotpr(workspace, 1, meanDifference, 1, &distance, dimensions);

    free(workspace);
    free(means);
    free(meanDifference);
    free(covarianceMatrix);
    free(covarianceMatrixInverse);

    return distance;
}

typedef struct ATSKalmanFilter
{
    double *Q;
    double *R;
    double *A;
    double *H;
    double *K;
    double *P;
    double *x;

    int stateDimensions;
    int observableDimensions;

} ATSKalmanFilter;

static void freeKalmanFilter(ATSKalmanFilter f)
{
    free(f.Q);
    free(f.R);
    free(f.A);
    free(f.H);
    free(f.K);
    free(f.P);
    free(f.x);
}

static ATSKalmanFilter ATSKalmanFilterMake(int stateDimensions, int observableDimensions)
{
    ATSKalmanFilter filter;

    filter.stateDimensions = stateDimensions;
    filter.observableDimensions = observableDimensions;

    filter.R = (double *)calloc(observableDimensions * observableDimensions, sizeof(double));
    filter.Q = (double *)calloc(stateDimensions * stateDimensions, sizeof(double));
    
    filter.A = (double *)calloc(stateDimensions * stateDimensions, sizeof(double));
    filter.H = (double *)calloc(stateDimensions * observableDimensions, sizeof(double));
    
    filter.K = (double *)calloc(stateDimensions * observableDimensions, sizeof(double));
    
    filter.P = (double *)calloc(stateDimensions * stateDimensions, sizeof(double));
    filter.x = (double *)calloc(stateDimensions, sizeof(double));

    return filter;
}

static ATSKalmanFilter ATSKalmanFilterMakeForGPS(double gpsError, double noise)
{
    ATSKalmanFilter filter = ATSKalmanFilterMake(4, 2);

    /** Setup the model */
    
    double d = 1.0;
    double A[16] = {
        1, 0, d, 0,
        0, 1, 0, d,
        0, 0, 1, 0,
        0, 0, 0, 1
    };

    double H[8] = {
        1, 0, 0, 0,
        0, 1, 0, 0
    };

    /** Setup the model error */

    double sig = gpsError * gpsError;
    double Q[16] = {
        sig, 0.0, 0.0, 0.0,
        0.0, sig, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0
    };

    double obs = sig * noise;
    double R[4] = {
        obs, 0.0,
        0.0, obs
    };

    /** Initial Guess */

    double inf = 1E12;
    double P[16] = {
        inf, 0.0, 0.0, 0.0,
        0.0, inf, 0.0, 0.0,
        0.0, 0.0, inf, 0.0,
        0.0, 0.0, 0.0, inf
    };

    vDSP_vaddD(&A[0], 1, filter.A, 1, filter.A, 1, 16);
    vDSP_vaddD(&H[0], 1, filter.H, 1, filter.H, 1, 8);
    vDSP_vaddD(&Q[0], 1, filter.Q, 1, filter.Q, 1, 16);
    vDSP_vaddD(&R[0], 1, filter.R, 1, filter.R, 1, 4);
    vDSP_vaddD(&P[0], 1, filter.P, 1, filter.P, 1, 16);

    return filter;
}

static void ATSKalmanFilterUpdate(ATSKalmanFilter *f, double *z)
{
    int sd = f->stateDimensions;
    int od = f->observableDimensions;

    /* Time update */

    // x[k] = A x[k-1]
    vDSP_mmulD(f->A, 1, f->x, 1, f->x, 1, sd, 1, sd);

    // AP = A P[k-1]
    double *AP = malloc(sd * sd * sizeof(double));
    vDSP_mmulD(f->A, 1, f->P, 1, AP, 1, sd, sd, sd);

    // AT = Transpose(A)
    double *AT = malloc(sd * sd * sizeof(double));
    vDSP_mtransD(f->A, 1, AT, 1, sd, sd);

    // D = A P[k-1] A^T = AP AT
    double *D = malloc(sd * sd * sizeof(double));
    vDSP_mmulD(AP, 1, AT, 1, D, 1, sd, sd, sd);

    // P[k] = A P[k-1] A^T + Q = D + Q
    vDSP_vaddD(D, 1, f->Q, 1, f->P, 1, sd * sd);

    free(AP);
    free(AT);
    free(D);

    /* Measurement Update */

    // Compute Gain K

    // HT = Transpose(H)
    double *HT = malloc(sd * od * sizeof(double));
    vDSP_mtransD(f->H, 1, HT, 1, sd, od);

    // PHT = P H^T = P HT
    double *PHT = malloc(sd * od * sizeof(double));
    vDSP_mmulD(f->P, 1, HT, 1, PHT, 1, sd, od, sd);

    // HPHT = H P H^T = H PHT
    double *HPHT = malloc(od * od * sizeof(double));
    vDSP_mmulD(f->H, 1, PHT, 1, HPHT, 1, od, od, sd);

    // G = H P H^T + R = HPHT + R
    double *G = malloc(od * od * sizeof(double));
    vDSP_vaddD(HPHT, 1, f->R, 1, G, 1, od * od);

    // IG = G^-1
    double *IG = ATSMatrixInverseD(G, od);

    // K = P H^T (H P H^T + R)^-1 = PHT * IG
    vDSP_mmulD(PHT, 1, IG, 1, f->K, 1, sd, od, od);

    free(HT);
    free(PHT);
    free(HPHT);
    free(G);
    free(IG);

    // Update x and P

    // Hx = H x
    double *Hx = malloc(od * sizeof(double));
    vDSP_mmulD(f->H, 1, f->x, 1, Hx, 1, od, 1, sd);

    // dz = z - Hx
    double *dz = malloc(od * sizeof(double));
    vDSP_vsubD(Hx, 1, z, 1, dz, 1, od);

    // Kdz = K dz
    double *Kdz = malloc(sd * sizeof(double));
    vDSP_mmulD(f->K, 1, dz, 1, Kdz, 1, sd, 1, od);

    // HP = H P
    double *HP = malloc(sd * od * sizeof(double));
    vDSP_mmulD(f->H, 1, f->P, 1, HP, 1, od, sd, sd);

    // KHP = K HP
    double *KHP = malloc(sd * sd * sizeof(double));
    vDSP_mmulD(f->K, 1, HP, 1, KHP, 1, sd, sd, od);

    // Update x
    // x = x + K dz
    vDSP_vaddD(f->x, 1, Kdz, 1, f->x, 1, sd);

    // Update P
    // P = P - KHP
    vDSP_vsubD(KHP, 1, f->P, 1, f->P, 1, sd * sd);

    free(Hx);
    free(dz);
    free(Kdz);
    free(HP);
    free(KHP);
}

static double *ATSKalmanFilterPredict(ATSKalmanFilter f)
{
    int od = f.observableDimensions;
    int sd = f.stateDimensions;
    
    double *Hx = malloc(od * sizeof(double));
    vDSP_mmulD(f.H, 1, f.x, 1, Hx, 1, od, 1, sd);
    
    return Hx;
}
