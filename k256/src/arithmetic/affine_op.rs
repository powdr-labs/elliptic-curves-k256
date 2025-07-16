use crate::AffinePoint;
use crate::FieldElement;
use core::ops::Add;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FieldBytes;
    use hex_literal::hex;
    #[test]
    fn test_addition() {
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

        let result = point1 + point2;
        assert_eq!(result.as_ref().unwrap().x, x3);
        assert_eq!(result.unwrap().y, y3);
    }
}
