export enum Sex { Male, Female }

export interface DosingStep {
    time: number; // seconds
    rate: number; // mcg/kg/min
    ce: number; // effect-site concentration (mcg/mL)
    cumulativeDoseMg: number; // mg
    phase: 'loading' | 'maintenance';
}

interface Params {
    drug: 'propofol' | 'remifentanil';
    concentrationMgPerMl: number;
    ageYears: number;
    sex: Sex;
    weightKg: number;
    heightCm: number;
    targetCe: number; // mcg/mL for propofol, mcg/mL for remifentanil
    intervalSec?: number; // simulation interval
    volumeMl: number; // infusion volume
    currentCe?: number; // starting effect-site concentration
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
        // assume concomitant opioid as in original repository logic
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
        // remifentanil (Eleveld-Remifentanil)
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

export function eleveldProtocol(params: Params): DosingStep[] {
    const intervalSec = params.intervalSec ?? 10;
    const startCe = params.currentCe ?? 0;
    const p = eleveldParameters(params);

    const steps: DosingStep[] = [];

    // initial bolus to reach targetCe if needed
    const deltaCe = Math.max(params.targetCe - startCe, 0); // mg/L == mcg/mL
    const bolusMg = deltaCe * p.vc;
    if (bolusMg > 0) {
        const bolusRate = (bolusMg / params.weightKg) * 60 * 1000; // mcg/kg/min
        steps.push({
            time: 0,
            rate: bolusRate,
            ce: startCe,
            cumulativeDoseMg: bolusMg,
            phase: 'loading'
        });
    }

    // infusion parameters
    const infusionRateMgMin = params.targetCe * p.cl1; // mg/min
    const ratePerKg = (infusionRateMgMin / params.weightKg) * 1000; // mcg/kg/min
    const totalInfDoseMg = params.volumeMl * params.concentrationMgPerMl;

    let C1 = startCe + deltaCe;
    let C2 = 0, C3 = 0, Ce = startCe + deltaCe;
    let cumulative = bolusMg;
    let infused = 0;

    const { k10, k12, k13, k21, k31, ke0 } = p;

    for (let t = intervalSec; ; t += intervalSec) {
        const dtMin = intervalSec / 60;
        const delivering = infused < totalInfDoseMg;
        const rateMgMin = delivering ? infusionRateMgMin : 0;

        const dC1 = (rateMgMin / p.vc) - (k10 + k12 + k13 + ke0) * C1 + k21 * C2 + k31 * C3;
        const dC2 = k12 * C1 - k21 * C2;
        const dC3 = k13 * C1 - k31 * C3;
        const dCe = ke0 * (C1 - Ce);

        C1 += dC1 * dtMin;
        C2 += dC2 * dtMin;
        C3 += dC3 * dtMin;
        Ce += dCe * dtMin;

        const added = rateMgMin * dtMin;
        infused += added;
        cumulative += added;

        steps.push({
            time: t,
            rate: delivering ? ratePerKg : 0,
            ce: Ce,
            cumulativeDoseMg: cumulative,
            phase: 'maintenance'
        });

        if (!delivering && Ce <= 0.05) break;
    }

    return steps;
}
