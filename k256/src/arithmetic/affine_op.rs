use crate::arithmetic::mul::{G1, G2, MINUS_B1, MINUS_B2, MINUS_LAMBDA, Radix16Decomposition};
use crate::arithmetic::projective::ENDOMORPHISM_BETA;
use crate::arithmetic::scalar::{Scalar, WideScalar};
use crate::{AffinePoint, FieldElement};
use core::ops::{Add, Mul, Neg};
use core::usize;
use elliptic_curve::group::prime::PrimeCurveAffine;
use elliptic_curve::point::AffineCoordinates;
use elliptic_curve::scalar::IsHigh;
use elliptic_curve::subtle::{Choice, ConditionallySelectable};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
/// Represents an ECC point
pub struct PowdrAffinePoint {
    pub x: FieldElement,
    pub y: FieldElement,
    pub infinity: u8,
}

impl From<AffinePoint> for PowdrAffinePoint {
    fn from(point: AffinePoint) -> Self {
        let x_bytes = point.x();
        let y_bytes = point.y();

        PowdrAffinePoint {
            x: FieldElement::from_bytes(&x_bytes).unwrap(),
            y: FieldElement::from_bytes(&y_bytes).unwrap(),
            infinity: if point.is_identity().into() { 1 } else { 0 },
        }
    }
}

impl Add<PowdrAffinePoint> for PowdrAffinePoint {
    type Output = PowdrAffinePoint;

    fn add(self, other: PowdrAffinePoint) -> PowdrAffinePoint {
        if self.infinity != 0 {
            return other;
        }
        if other.infinity != 0 {
            return self;
        }

        if other.x == self.x {
            if self.y == other.y {
                return self.double();
            } else {
                //  x1 == x2 but y1 != y2 → vertical line → point at infinity
                return PowdrAffinePoint::IDENTITY;
            }
        }

        // normalization is needed here to ensure the initial magnitude is 1.
        let x1 = self.x.normalize();
        let x2 = other.x.normalize();
        let y1 = self.y.normalize();
        let y2 = other.y.normalize();

        let dx = (x2 - x1).normalize();
        let invert = dx.invert().unwrap();

        let dy = y2 - y1;
        let lambda = dy * invert;

        assert_eq!(
            FieldElement::from_u64(1).normalize(),
            (invert * dx).normalize()
        );

        let x3 = lambda.square() - x1 - x2;
        let y3 = lambda * (x1 + x3.negate(5)) - y1;

        PowdrAffinePoint {
            x: x3.normalize(),
            y: y3.normalize(),
            infinity: 0,
        }
    }
}

impl Neg for PowdrAffinePoint {
    type Output = PowdrAffinePoint;

    fn neg(self) -> PowdrAffinePoint {
        PowdrAffinePoint::neg(&self)
    }
}

impl Mul<Scalar> for PowdrAffinePoint {
    type Output = PowdrAffinePoint;

    fn mul(self, other: Scalar) -> PowdrAffinePoint {
        mul(&self, &other)
    }
}

impl PowdrAffinePoint {
    /// Additive identity of the group: the point at infinity.
    pub const IDENTITY: Self = Self {
        x: FieldElement::ZERO,
        y: FieldElement::ZERO,
        infinity: 1,
    };

    /// Base point of secp256k1.
    pub const GENERATOR: Self = Self {
        x: AffinePoint::GENERATOR.x,
        y: AffinePoint::GENERATOR.y,
        infinity: 0,
    };

    /// Double the point.
    pub fn double(self) -> PowdrAffinePoint {
        let x = self.x.normalize();
        let y = self.y.normalize();

        if y.is_zero().into() {
            return PowdrAffinePoint::IDENTITY;
        }

        let num = FieldElement::from(3u64) * x.square();
        let denom = (FieldElement::from(2u64) * y).normalize();
        let lambda = num * denom.invert().unwrap();

        let x3 = lambda.square() - FieldElement::from(2u64) * x;
        let y3 = lambda * (x + x3.negate(3)) - y;

        PowdrAffinePoint {
            x: x3.normalize(),
            y: y3.normalize(),
            infinity: 0,
        }
    }

    fn neg(&self) -> Self {
        PowdrAffinePoint {
            x: self.x,
            y: self.y.negate(1).normalize_weak(),
            infinity: self.infinity,
        }
    }

    /// Calculates `k * G`, where `G` is the generator, using precomputed tables.
    // TODO: can use precomputed tables for better performance
    pub(super) fn mul_by_generator(k: &Scalar) -> PowdrAffinePoint {
        PowdrAffinePoint::GENERATOR * *k
    }

    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        PowdrAffinePoint {
            x: FieldElement::conditional_select(&a.x, &b.x, choice),
            y: FieldElement::conditional_select(&a.y, &b.y, choice),
            infinity: choice.unwrap_u8(),
        }
    }

    /// Calculates SECP256k1 endomorphism: `self * lambda`.
    pub fn endomorphism(&self) -> Self {
        Self {
            x: self.x * &ENDOMORPHISM_BETA,
            y: self.y,
            infinity: self.infinity,
        }
    }
}

#[inline(always)]
fn mul(x: &PowdrAffinePoint, k: &Scalar) -> PowdrAffinePoint {
    lincomb(&[(*x, *k)])
}

pub fn lincomb<const N: usize>(
    points_and_scalars: &[(PowdrAffinePoint, Scalar); N],
) -> PowdrAffinePoint {
    let mut tables = [(LookupTable::default(), LookupTable::default()); N];
    let mut digits = [(
        Radix16Decomposition::<33>::default(),
        Radix16Decomposition::<33>::default(),
    ); N];

    lincomb_pippenger(points_and_scalars, &mut tables, &mut digits)
}

fn lincomb_pippenger(
    xks: &[(PowdrAffinePoint, Scalar)],
    tables: &mut [(LookupTable, LookupTable)],
    digits: &mut [(Radix16Decomposition<33>, Radix16Decomposition<33>)],
) -> PowdrAffinePoint {
    xks.iter().enumerate().for_each(|(i, (x, k))| {
        let (r1, r2) = decompose_scalar(k);
        let x_beta = x.endomorphism();
        let (r1_sign, r2_sign) = (r1.is_high(), r2.is_high());

        let (r1_c, r2_c) = (
            Scalar::conditional_select(&r1, &-r1, r1_sign),
            Scalar::conditional_select(&r2, &-r2, r2_sign),
        );

        tables[i] = (
            LookupTable::from(&PowdrAffinePoint::conditional_select(x, &-*x, r1_sign)),
            LookupTable::from(&PowdrAffinePoint::conditional_select(
                &x_beta, &-x_beta, r2_sign,
            )),
        );

        digits[i] = (
            Radix16Decomposition::<33>::new(&r1_c),
            Radix16Decomposition::<33>::new(&r2_c),
        )
    });

    let mut acc = PowdrAffinePoint::IDENTITY;
    for component in 0..xks.len() {
        let (digit1, digit2) = digits[component];
        let (table1, table2) = tables[component];

        acc = table1.select(digit1.0[32]) + acc;
        acc = table2.select(digit2.0[32]) + acc;
    }

    for i in (0..32).rev() {
        for _j in 0..4 {
            acc = acc.double();
        }

        for component in 0..xks.len() {
            let (digit1, digit2) = digits[component];
            let (table1, table2) = tables[component];

            acc = table1.select(digit1.0[i]) + acc;
            acc = table2.select(digit2.0[i]) + acc;
        }
    }

    acc
}

/// Find r1 and r2 given k, such that r1 + r2 * lambda == k mod n.
fn decompose_scalar(k: &Scalar) -> (Scalar, Scalar) {
    // these _vartime calls are constant time since the shift amount is constant
    let c1 = WideScalar::mul_shift_vartime(k, &G1, 384) * MINUS_B1;
    let c2 = WideScalar::mul_shift_vartime(k, &G2, 384) * MINUS_B2;
    let r2 = c1 + c2;
    let r1 = k + r2 * MINUS_LAMBDA;

    (r1, r2)
}

#[derive(Copy, Clone, Default)]
struct LookupTable([PowdrAffinePoint; 8]);

impl From<&PowdrAffinePoint> for LookupTable {
    fn from(p: &PowdrAffinePoint) -> Self {
        let mut points = [*p; 8];
        for j in 0..7 {
            points[j + 1] = p.clone() + points[j].clone();
        }
        LookupTable(points)
    }
}

impl LookupTable {
    /// Given -8 <= x <= 8, returns x * p in constant time.
    fn select(&self, x: i8) -> PowdrAffinePoint {
        debug_assert!((-8..=8).contains(&x));

        if x == 0 {
            PowdrAffinePoint::IDENTITY
        } else {
            let abs = x.abs() as usize;
            let mut point = self.0[abs - 1];

            if x < 0 {
                point.y = -point.y;
            }

            point
        }
    }
}

#[cfg(test)]
mod tests {
    use std::println;

    use super::*;
    use crate::FieldBytes;
    use crate::arithmetic::{ProjectivePoint, Scalar};
    use elliptic_curve::PrimeField;
    use elliptic_curve::{
        Field, Group,
        rand_core::{OsRng, TryRngCore},
    };
    use hex_literal::hex;

    #[test]
    fn test_addition_double() {
        let x1: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[{
                let mut bytes = [0u8; 32];
                bytes[31] = 1;
                bytes
            }])[0],
        )
        .unwrap();
        let y1: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[<[_; 32]>::try_from(hex!(
                "4218F20AE6C646B363DB68605822FB14264CA8D2587FDD6FBC750D587E76A7EE"
            ))
            .unwrap()])[0],
        )
        .unwrap()
        .normalize();

        let x2: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[{
                let mut bytes = [0u8; 32];
                bytes[31] = 2;
                bytes
            }])[0],
        )
        .unwrap();

        let y2: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[<[_; 32]>::try_from(hex!(
                "990418D84D45F61F60A56728F5A10317BDB3A05BDA4425E3AEE079F8A847A8D1"
            ))
            .unwrap()])[0],
        )
        .unwrap()
        .normalize();

        let x3: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[<[_; 32]>::try_from(hex!(
                "F23A2D865C24C99CC9E7B99BD907FB93EBD6CCCE106BCCCB0082ACF8315E67BE"
            ))
            .unwrap()])[0],
        )
        .unwrap()
        .normalize();

        let y3: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[<[_; 32]>::try_from(hex!(
                "791DFC78B49C9B5882867776F18BA7883ED0BAE1C0A856D26D41D38FB47345B4"
            ))
            .unwrap()])[0],
        )
        .unwrap()
        .normalize();

        let x4: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[<[_; 32]>::try_from(hex!(
                "33333333333333333333333333333333333333333333333333333332FFFFFF3B"
            ))
            .unwrap()])[0],
        )
        .unwrap()
        .normalize();

        let y4: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[<[_; 32]>::try_from(hex!(
                "3916485F2C3D80C62048C6FD8ACBF71EED11987A55CC10ABDC4E4A25C4EC54AC"
            ))
            .unwrap()])[0],
        )
        .unwrap()
        .normalize();

        let point1 = PowdrAffinePoint {
            x: x1,
            y: y1,
            infinity: 0,
        };
        let point2 = PowdrAffinePoint {
            x: x2,
            y: y2,
            infinity: 0,
        };

        let addition = point1 + point2.clone();
        let double = point2.double();
        assert_eq!(addition.x, x3);
        assert_eq!(addition.y, y3);

        assert_eq!(double.x, x4.normalize());
        assert_eq!(double.y, y4.normalize());
    }

    #[test]
    fn test_multiplication() {
        let x1: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[{
                let mut bytes = [0u8; 32];
                bytes[31] = 1;
                bytes
            }])[0],
        )
        .unwrap();
        let y1: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[<[_; 32]>::try_from(hex!(
                "4218F20AE6C646B363DB68605822FB14264CA8D2587FDD6FBC750D587E76A7EE"
            ))
            .unwrap()])[0],
        )
        .unwrap()
        .normalize();

        let x5: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[<[_; 32]>::try_from(hex!(
                "6D6D216817A448DC312FEE586FA306D189CB404A9CAF72D90308797F38934A19"
            ))
            .unwrap()])[0],
        )
        .unwrap()
        .normalize();

        let y5: FieldElement = FieldElement::from_bytes(
            &FieldBytes::cast_slice_from_core(&[<[_; 32]>::try_from(hex!(
                "2C9BB19372B2E1B830B5F4D92ADBAFEAAEB612026122E571D1BEA76D742F279E"
            ))
            .unwrap()])[0],
        )
        .unwrap()
        .normalize();

        let scalar = Scalar::from_u128(12345678);

        let point1 = PowdrAffinePoint {
            x: x1,
            y: y1,
            infinity: 0,
        };
        let point5 = PowdrAffinePoint {
            x: x5,
            y: y5,
            infinity: 0,
        };

        let multiplication = point1 * scalar;
        assert_eq!(multiplication.x, point5.x);
        assert_eq!(multiplication.y, point5.y);
    }

    #[test]
    fn test_lincomb() {
        let x =
            PowdrAffinePoint::from(ProjectivePoint::random(&mut OsRng.unwrap_mut()).to_affine());
        let y =
            PowdrAffinePoint::from(ProjectivePoint::random(&mut OsRng.unwrap_mut()).to_affine());
        let k = Scalar::random(&mut OsRng.unwrap_mut());
        let l = Scalar::random(&mut OsRng.unwrap_mut());

        println!("k: {:?}", k);
        println!("x: {:?}", x);
        println!("y: {:?}", y);
        println!("l: {:?}", l);

        let reference = x * k + y * l;
        let test = lincomb(&[(x, k), (y, l)]);
        assert_eq!(reference, test);
    }
}
