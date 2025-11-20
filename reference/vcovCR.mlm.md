# Cluster-robust variance-covariance matrix for an mlm object.

`vcovCR` returns a sandwich estimate of the variance-covariance matrix
of a set of regression coefficient estimates from an `mlm` object.

## Usage

``` r
# S3 method for class 'mlm'
vcovCR(obj, cluster, type, target, inverse_var, form = "sandwich", ...)
```

## Arguments

- obj:

  Fitted model for which to calculate the variance-covariance matrix

- cluster:

  Optional expression or vector indicating which observations belong to
  the same cluster. If not specified, each row of the data will be
  treated as a separate cluster.

- type:

  Character string specifying which small-sample adjustment should be
  used, with available options `"CR0"`, `"CR1"`, `"CR1p"`, `"CR1S"`,
  `"CR2"`, or `"CR3"`. See "Details" section of
  [`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
  for further information.

- target:

  Optional matrix or vector describing the working variance-covariance
  model used to calculate the `CR2` and `CR4` adjustment matrices. If
  not specified, the target is taken to be an identity matrix.

- inverse_var:

  Optional logical indicating whether the weights used in fitting the
  model are inverse-variance. If not specified, `vcovCR` will attempt to
  infer a value.

- form:

  Controls the form of the returned matrix. The default `"sandwich"`
  will return the sandwich variance-covariance matrix. Alternately,
  setting `form = "meat"` will return only the meat of the sandwich and
  setting `form = B`, where `B` is a matrix of appropriate dimension,
  will return the sandwich variance-covariance matrix calculated using
  `B` as the bread. `form = "estfun"` will return the (appropriately
  scaled) estimating function, the transposed crossproduct of which is
  equal to the sandwich variance-covariance matrix.

- ...:

  Additional arguments available for some classes of objects.

## Value

An object of class `c("vcovCR","clubSandwich")`, which consists of a
matrix of the estimated variance of and covariances between the
regression coefficient estimates.

## See also

[`vcovCR`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)

## Examples

``` r
iris_fit <- lm(cbind(Sepal.Length, Sepal.Width) ~ Species + 
               Petal.Length + Petal.Width, data = iris)
Vcluster <- vcovCR(iris_fit, type = "CR2")
Vcluster
#>                                Sepal.Length:(Intercept)
#> Sepal.Length:(Intercept)                    0.012547708
#> Sepal.Length:Speciesversicolor              0.013463389
#> Sepal.Length:Speciesvirginica               0.022738726
#> Sepal.Length:Petal.Length                  -0.007227063
#> Sepal.Length:Petal.Width                    0.002810727
#> Sepal.Width:(Intercept)                     0.006324275
#> Sepal.Width:Speciesversicolor               0.004545412
#> Sepal.Width:Speciesvirginica                0.007262239
#> Sepal.Width:Petal.Length                   -0.003231810
#> Sepal.Width:Petal.Width                     0.002072658
#>                                Sepal.Length:Speciesversicolor
#> Sepal.Length:(Intercept)                         0.0134633889
#> Sepal.Length:Speciesversicolor                   0.0380390801
#> Sepal.Length:Speciesvirginica                    0.0579846604
#> Sepal.Length:Petal.Length                       -0.0087024291
#> Sepal.Length:Petal.Width                        -0.0109399201
#> Sepal.Width:(Intercept)                          0.0045454124
#> Sepal.Width:Speciesversicolor                    0.0125587855
#> Sepal.Width:Speciesvirginica                     0.0176278866
#> Sepal.Width:Petal.Length                        -0.0039842873
#> Sepal.Width:Petal.Width                         -0.0000394729
#>                                Sepal.Length:Speciesvirginica
#> Sepal.Length:(Intercept)                         0.022738726
#> Sepal.Length:Speciesversicolor                   0.057984660
#> Sepal.Length:Speciesvirginica                    0.096974290
#> Sepal.Length:Petal.Length                       -0.013179441
#> Sepal.Length:Petal.Width                        -0.021314939
#> Sepal.Width:(Intercept)                          0.007262239
#> Sepal.Width:Speciesversicolor                    0.017627887
#> Sepal.Width:Speciesvirginica                     0.026715559
#> Sepal.Width:Petal.Length                        -0.005518284
#> Sepal.Width:Petal.Width                         -0.001400422
#>                                Sepal.Length:Petal.Length
#> Sepal.Length:(Intercept)                    -0.007227063
#> Sepal.Length:Speciesversicolor              -0.008702429
#> Sepal.Length:Speciesvirginica               -0.013179441
#> Sepal.Length:Petal.Length                    0.005876684
#> Sepal.Length:Petal.Width                    -0.006412572
#> Sepal.Width:(Intercept)                     -0.003231810
#> Sepal.Width:Speciesversicolor               -0.003984287
#> Sepal.Width:Speciesvirginica                -0.005518284
#> Sepal.Width:Petal.Length                     0.002335906
#> Sepal.Width:Petal.Width                     -0.001991957
#>                                Sepal.Length:Petal.Width Sepal.Width:(Intercept)
#> Sepal.Length:(Intercept)                   0.0028107268             0.006324275
#> Sepal.Length:Speciesversicolor            -0.0109399201             0.004545412
#> Sepal.Length:Speciesvirginica             -0.0213149395             0.007262239
#> Sepal.Length:Petal.Length                 -0.0064125723            -0.003231810
#> Sepal.Length:Petal.Width                   0.0272498285             0.002072658
#> Sepal.Width:(Intercept)                    0.0020726578             0.013117244
#> Sepal.Width:Speciesversicolor             -0.0000394729             0.013714416
#> Sepal.Width:Speciesvirginica              -0.0014004219             0.019456928
#> Sepal.Width:Petal.Length                  -0.0019919567            -0.007453511
#> Sepal.Width:Petal.Width                    0.0049052390             0.003506130
#>                                Sepal.Width:Speciesversicolor
#> Sepal.Length:(Intercept)                        0.0045454124
#> Sepal.Length:Speciesversicolor                  0.0125587855
#> Sepal.Length:Speciesvirginica                   0.0176278866
#> Sepal.Length:Petal.Length                      -0.0039842873
#> Sepal.Length:Petal.Width                       -0.0000394729
#> Sepal.Width:(Intercept)                         0.0137144162
#> Sepal.Width:Speciesversicolor                   0.0349169656
#> Sepal.Width:Speciesvirginica                    0.0472678537
#> Sepal.Width:Petal.Length                       -0.0105186360
#> Sepal.Width:Petal.Width                        -0.0025487274
#>                                Sepal.Width:Speciesvirginica
#> Sepal.Length:(Intercept)                        0.007262239
#> Sepal.Length:Speciesversicolor                  0.017627887
#> Sepal.Length:Speciesvirginica                   0.026715559
#> Sepal.Length:Petal.Length                      -0.005518284
#> Sepal.Length:Petal.Width                       -0.001400422
#> Sepal.Width:(Intercept)                         0.019456928
#> Sepal.Width:Speciesversicolor                   0.047267854
#> Sepal.Width:Speciesvirginica                    0.069056876
#> Sepal.Width:Petal.Length                       -0.013526991
#> Sepal.Width:Petal.Width                        -0.007668001
#>                                Sepal.Width:Petal.Length Sepal.Width:Petal.Width
#> Sepal.Length:(Intercept)                   -0.003231810            0.0020726578
#> Sepal.Length:Speciesversicolor             -0.003984287           -0.0000394729
#> Sepal.Length:Speciesvirginica              -0.005518284           -0.0014004219
#> Sepal.Length:Petal.Length                   0.002335906           -0.0019919567
#> Sepal.Length:Petal.Width                   -0.001991957            0.0049052390
#> Sepal.Width:(Intercept)                    -0.007453511            0.0035061303
#> Sepal.Width:Speciesversicolor              -0.010518636           -0.0025487274
#> Sepal.Width:Speciesvirginica               -0.013526991           -0.0076680010
#> Sepal.Width:Petal.Length                    0.005576498           -0.0042578708
#> Sepal.Width:Petal.Width                    -0.004257871            0.0131707991
```
