use nalgebra::{self as na, dvector};
use na::{DMatrix, DVector};

fn mahalanobis_distance(x: &DVector<f64>, d: &DMatrix<f64>) -> f64 {

    let mu = mean_vector(d);

    // Calculer la différence entre le point x et la moyenne mu
    let diff = x - mu;

    let covariance_inv = covariance_matrix(d).try_inverse().expect("La matrice de covariance doit être inversible");

    // Calculer le produit (diff ^ T) * covariance_inv
    let product = covariance_inv * diff.clone();

    // Calculer la distance de Mahalanobis : sqrt(diff^T * covariance_inv * diff)
    let dist = diff.dot(&product).sqrt();
    dist
}

// Fonction pour calculer la moyenne (vecteur) d'une matrice
fn mean_vector(d: &DMatrix<f64>) -> DVector<f64> {

    d.row_mean_tr()
}

// Fonction pour calculer la matrice de covariance d'une distribution D
fn covariance_matrix(d: &DMatrix<f64>) -> DMatrix<f64> {
    let n = d.nrows() as f64;

    let mean = mean_vector(d);

    let mut covariance = DMatrix::zeros(d.ncols(), d.ncols());
    
    for i in 0..d.nrows() {
        
        let diff = &d.row(i).transpose() - mean.clone();
        covariance += diff.clone() * diff.transpose();
    }
    
    covariance / (n - 1.0)
}

pub fn main() {
    // Exemple de données
    let x = dvector![66.0, 640.0, 44.0];

    let d = DMatrix::from_row_slice(5,3, &[64.0, 580.0, 29.0, 
                                                                                        66.0, 570.0, 33.0,
                                                                                        68.0, 590.0, 37.0,
                                                                                        69.0, 660.0, 46.0,
                                                                                        73.0, 600.0, 55.0]);


    let covariance = covariance_matrix(&d);
    println!("Matrice de covariance :\n{}", covariance);

    // Calcul de la distance de Mahalanobis
    let distance = mahalanobis_distance(&x, &d);

    println!("La distance de Mahalanobis est : {}", distance);
}

