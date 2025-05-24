export enum Sex { Male, Female }

export interface DosingStep {
    time: number;        // seconds
    rate: number;        // mcg/kg/min
    ce: number;          // effect-site concentration (mcg/mL)
    cumulativeDoseMg: number; // mg
    phase: 'loading' | 'maintenance' | 'decay';
}

interface Params {
    drug: 'propofol' | 'remifentanil';
    concentrationMgPerMl: number;
    ageYears: number;
    sex: Sex;
    weightKg: number;
    heightCm: number;
    targetCe: number;
    intervalSec?: number;
    volumeMl?: number;          // total infusion volume after loading
    currentCe?: number;         // starting effect-site concentration
    loadingDoseMg?: number;     // optional precomputed loading dose
    loadingDurationSec?: number;// optional loading infusion duration
}

function fsigmoid(x: number, y: number, z: number): number {
    return Math.pow(x, z) / (Math.pow(x, z) + Math.pow(y, z));
}

function fageing(age: number, x: number): number {
    return Math.exp(x * (age - 35));
}

function fcentral(x: number): number {
    return fsigmoid(x, 33.6, 1);
}

function fclmaturation(weeks: number): number {
    return fsigmoid(weeks, 42.3, 9.06);
}

function fq3maturation(weeks: number): number {
    return fsigmoid(weeks + 40, 68.3, 1);
}

function fffm(weight: number, height: number, age: number, sex: Sex): number {
    const bmi = weight / Math.pow(height / 100, 2);
    if (sex === Sex.Male) {
        return (0.88 + (1 - 0.88) / (1 + Math.pow(age / 13.4, -12.7))) * ((9270 * weight) / (6680 + 216 * bmi));
    }
    return (1.11 + (1 - 1.11) / (1 + Math.pow(age / 7.1, -1.1))) * ((9270 * weight) / (8780 + 244 * bmi));
}

function eleveldParameters(p: Params) {
    const toWeeks = 52.1429;
    let k10: number, k12: number, k13: number, k21: number, k31: number, ke0: number;
    let vc = 0, cl1 = 0, cl2 = 0, cl3 = 0, v2 = 0, v3 = 0;

    if (p.drug === 'propofol') {
        const PMA = p.ageYears * toWeeks + 40;
        vc = 6.28 * fcentral(p.weightKg) / fcentral(70);
        v2 = 25.5 * p.weightKg / 70 * fageing(p.ageYears, -0.0156);
        const ffmRef = (0.88 + (1 - 0.88) / (1 + Math.pow(35 / 13.4, -12.7))) * ((9270 * 70) / (6680 + 216 * 24.22145));
        v3 = 273 * fffm(p.weightKg, p.heightCm, p.ageYears, p.sex) / ffmRef;
        v3 *= Math.exp(-0.0138 * p.ageYears);
        const v2ref = 25.5;
        const v3ref = 273;
        if (p.sex === Sex.Male) {
            cl1 = 1.79 * Math.pow(p.weightKg / 70, 0.75) * (fclmaturation(PMA) / fclmaturation(35 * toWeeks + 40));
        } else {
            cl1 = 2.1 * Math.pow(p.weightKg / 70, 0.75) * (fclmaturation(PMA) / fclmaturation(35 * toWeeks + 40));
        }
        cl2 = 1.75 * Math.pow(v2 / v2ref, 0.75) * (1 + 1.3 * (1 - fq3maturation(p.ageYears * toWeeks)));
        cl3 = 1.11 * Math.pow(v3 / v3ref, 0.75) * (fq3maturation(p.ageYears * toWeeks) / fq3maturation(35 * toWeeks));
        ke0 = 0.146 * Math.pow(p.weightKg / 70, -0.25);
        k10 = cl1 / vc; k12 = cl2 / vc; k13 = cl3 / vc; k21 = cl2 / v2; k31 = cl3 / v3;
    } else {
        const ffmRef = (0.88 + (1 - 0.88) / (1 + Math.pow(35 / 13.4, -12.7))) * ((9270 * 70) / (6680 + 216 * 24.22145));
        const size = fffm(p.weightKg, p.heightCm, p.ageYears, p.sex) / ffmRef;
        const ksex = p.sex === Sex.Male ? 1 : 1 + 0.47 * fsigmoid(p.ageYears, 12, 6) * (1 - fsigmoid(p.ageYears, 45, 6));
        vc = 5.81 * size * fageing(p.ageYears, -0.00554);
        v2 = 8.82 * size * fageing(p.ageYears, -0.00327) * ksex;
        v3 = 5.03 * size * fageing(p.ageYears, -0.0315) * Math.exp(-0.026 * (p.weightKg - 70));
        cl1 = 2.58 * Math.pow(size, 0.75) * (fsigmoid(p.weightKg, 2.88, 2) / fsigmoid(70, 2.88, 2)) * ksex * fageing(p.ageYears, -0.00327);
        cl2 = 1.72 * Math.pow(v2 / 8.82, 0.75) * fageing(p.ageYears, -0.00554) * ksex;
        cl3 = 0.124 * Math.pow(v3 / 5.03, 0.75) * fageing(p.ageYears, -0.00554);
        ke0 = p.ageYears <= 16 ? 0.71 : 1.09 * fageing(p.ageYears, -0.0289);
        k10 = cl1 / vc; k12 = cl2 / vc; k13 = cl3 / vc; k21 = cl2 / v2; k31 = cl3 / v3;
    }
    return { vc, v2, v3, k10, k12, k13, k21, k31, ke0, cl1 };
}

export function calculateLoadingDose(params: Params) {
    const p = eleveldParameters(params);
    const startCe = params.currentCe ?? 0;
    const deltaCe = Math.max(params.targetCe - startCe, 0);
    const doseMg = deltaCe * p.vc;
    const mgPerSec = params.concentrationMgPerMl * (1000 / 3600);
    const durationSec = doseMg / mgPerSec;
    return { loadingDoseMg: doseMg, durationSec };
}

export function eleveldProtocol(params: Params): DosingStep[] {
    const interval = params.intervalSec ?? 10;
    const p = eleveldParameters(params);

    const maxRateMgSec = params.concentrationMgPerMl * (1000 / 3600);

    let { loadingDoseMg, loadingDurationSec } = params;
    if (loadingDoseMg === undefined || loadingDurationSec === undefined) {
        const calc = calculateLoadingDose(params);
        if (loadingDoseMg === undefined) loadingDoseMg = calc.loadingDoseMg;
        if (loadingDurationSec === undefined) loadingDurationSec = calc.durationSec;
    }

    // ensure loading rate does not exceed pump capability
    let bolusRateMgSec = loadingDoseMg / loadingDurationSec;
    if (bolusRateMgSec > maxRateMgSec) {
        bolusRateMgSec = maxRateMgSec;
        loadingDurationSec = loadingDoseMg / bolusRateMgSec;
    }

    const steps: DosingStep[] = [];

    let C1 = params.currentCe ?? 0;
    let C2 = 0, C3 = 0, Ce = params.currentCe ?? 0;
    let cumulative = 0;
    let time = 0;

    const { k10, k12, k13, k21, k31, ke0 } = p;

    // integrate loading infusion
    while (time < loadingDurationSec) {
        const dt = Math.min(interval, loadingDurationSec - time);
        const ratePerKg = (bolusRateMgSec * 60) / params.weightKg * 1000;

        const dC1 = (bolusRateMgSec / p.vc) - (k10 + k12 + k13) * C1 + k21 * C2 + k31 * C3;
        const dC2 = k12 * C1 - k21 * C2;
        const dC3 = k13 * C1 - k31 * C3;
        const dCe = ke0 * (C1 - Ce);

        C1 += dC1 * dt;
        C2 += dC2 * dt;
        C3 += dC3 * dt;
        Ce += dCe * dt;

        cumulative += bolusRateMgSec * dt;
        time += dt;

        steps.push({
            time: Math.round(time),
            rate: ratePerKg,
            ce: Ce,
            cumulativeDoseMg: cumulative,
            phase: 'loading'
        });
    }

    let remainingMg = (params.volumeMl ?? 0) * params.concentrationMgPerMl;
    while (remainingMg > 0 || Ce > 0.01) {
        const dt = interval;
        let infusionRateMgSec = 0;
        if (remainingMg > 0) {
            const error = params.targetCe - Ce;
            let desiredMgSec = (params.targetCe * p.cl1 + error * p.vc * ke0) / 60;
            if (desiredMgSec < 0) desiredMgSec = 0;
            infusionRateMgSec = Math.min(maxRateMgSec, desiredMgSec);
            const maxDeliver = remainingMg / dt;
            if (infusionRateMgSec > maxDeliver) {
                infusionRateMgSec = maxDeliver;
            }
            remainingMg -= infusionRateMgSec * dt;
        }

        const ratePerKg = (infusionRateMgSec * 60) / params.weightKg * 1000;

        const dC1 = (infusionRateMgSec / p.vc) - (k10 + k12 + k13) * C1 + k21 * C2 + k31 * C3;
        const dC2 = k12 * C1 - k21 * C2;
        const dC3 = k13 * C1 - k31 * C3;
        const dCe = ke0 * (C1 - Ce);

        C1 += dC1 * dt;
        C2 += dC2 * dt;
        C3 += dC3 * dt;
        Ce += dCe * dt;

        cumulative += infusionRateMgSec * dt;
        time += dt;

        steps.push({
            time: Math.round(time),
            rate: ratePerKg,
            ce: Ce,
            cumulativeDoseMg: cumulative,
            phase: remainingMg > 0 ? 'maintenance' : 'decay'
        });

        if (remainingMg <= 0 && Ce <= 0.01) break;
    }

    return steps;
}
