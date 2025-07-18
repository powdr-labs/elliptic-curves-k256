use crate::FieldElement;
use crate::arithmetic::scalar::{Scalar, WideScalar};
use core::ops::Add;
use elliptic_curve::{ops::LinearCombination, scalar::IsHigh};

#[derive(Clone, Copy, Debug, PartialEq, Eq,Default)]
pub struct PowdrAffinePoint {
    x: FieldElement,
    y: FieldElement,
    infinity: bool,
}

impl Add<PowdrAffinePoint> for PowdrAffinePoint {
    type Output = Option<PowdrAffinePoint>;

    fn add(self, other: PowdrAffinePoint) -> Option<PowdrAffinePoint> {
        self.x.normalize();
        self.y.normalize();
        other.x.normalize();
        other.y.normalize();

        let dx = (other.x - self.x).normalize();
        if dx.is_zero().into() {
            return None;
        }
        let invert = dx.invert().unwrap();

        let dy = other.y - self.y;
        let lambda = dy * invert;

        assert_eq!(
            FieldElement::from_u64(1).normalize(),
            (invert * dx).normalize()
        );

        let x3 = lambda.square() - self.x - other.x;
        let y3 = lambda * (self.x + x3.negate(5)) - self.y;

        Some(PowdrAffinePoint {
            x: x3.normalize(),
            y: y3.normalize(),
            infinity: false,
        })
    }
}

impl PowdrAffinePoint {
    pub fn double(self) -> Option<PowdrAffinePoint> {
        let x = self.x.normalize();
        let y = self.y.normalize();

        if y.is_zero().into() {
            return None;
        }

        let num = FieldElement::from(3u64) * x.square();
        let denom = (FieldElement::from(2u64) * y).normalize();
        let lambda = num * denom.invert().unwrap();

        let x3 = lambda.square() - FieldElement::from(2u64) * x;
        let y3 = lambda * (x + x3.negate(3)) - y;

        Some(PowdrAffinePoint {
            x: x3.normalize(),
            y: y3.normalize(),
            infinity: false,
        })
    }
}

fn lincomb(points_and_scalars: &[(PowdrAffinePoint, Scalar); N]) -> PowdrAffinePoint {
    let mut tables = [(LookupTable::default(), LookupTable::default()); N];
    let mut digits = [(
        Radix16Decomposition::<33>::default(),
        Radix16Decomposition::<33>::default(),
    ); N];

    lincomb(points_and_scalars, &mut tables, &mut digits)
}

#[derive(Copy, Clone, Default)]
struct LookupTable([PowdrAffinePoint; 8]);

impl From<&PowdrAffinePoint> for LookupTable {
    fn from(p: &PowdrAffinePoint) -> Self {
        let mut points = [*p; 8];
        for j in 0..7 {
            points[j + 1] = (p.clone() + points[j].clone()).unwrap();
        }
        LookupTable(points)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FieldBytes;
    use hex_literal::hex;
    use sha2::digest::typenum::Double;
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
            infinity: false,
        };
        let point2 = PowdrAffinePoint {
            x: x2,
            y: y2,
            infinity: false,
        };

        let addition = point1 + point2.clone();
        let double = point2.double();
        assert_eq!(addition.as_ref().unwrap().x, x3);
        assert_eq!(addition.unwrap().y, y3);

        assert_eq!(double.as_ref().unwrap().x, x4.normalize());
        assert_eq!(double.unwrap().y, y4.normalize());
    }
}
