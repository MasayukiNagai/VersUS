<!DOCTYPE html>
<html lang="en">

<head>
  <!-- Required meta tags -->
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="An interface for exploring variants of uncertain significance">
  <meta name="author" content="Masayuki Nagai">

  <link href="css/bootstrap.min.css" rel="stylesheet">

  <title>VersUS(beta)</title>

  <style type="text/css">
    /* body { padding-top: 70px; } */
  </style>

  <!-- <style type="text/css">
    .navbar-custom {
    background-color: #633974 ;
    }
    /* change the brand and text color */
    .navbar-custom .navbar-brand,
    .navbar-custom .navbar-text {
    color: rgba(255,255,255,.8);
    }

    /* change the color of active or hovered links */
    .navbar-custom .nav-item.active .nav-link,
    .navbar-custom .nav-item:hover .nav-link {
    color: #ffffff;
    }
  </style> -->
</head>

<body>
<nav class="navbar navbar-default navbar-expand-md navbar-light bg-light navbar-fixed-top">
  <div class="container">
    <a class="navbar-brand" href="#">VersUS</a>
    <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarSupportedContent"
      aria-controls="navbarSupportedContent" aria-expanded="false" aria-label="Toggle navigation">
      <span class="navbar-toggler-icon"></span>
    </button>

    <div class="collapse navbar-collapse" id="navbarSupportedContent">
      <ul class="navbar-nav me-auto">
        <li class="nav-item active">
          <a class="nav-link" href="#">Home</a>
        </li>
        <li class="nav-item">
          <a class="nav-link" href="#">Tree view</a>
        </li>
        <!-- <li class="nav-item dropdown">
          <a class="nav-link dropdown-toggle" href="#" id="navbarDropdownMenuLink" role="button" data-bs-toggle="dropdown"
            aria-expanded="false">
            Dropdown link
          </a>
          <ul class="dropdown-menu" aria-labelledby="navbarDropdownMenuLink">
            <li><a class="dropdown-item" href="#">Action</a></li>
            <li><a class="dropdown-item" href="#">Another action</a></li>
            <li><a class="dropdown-item" href="#">Something else here</a></li>
          </ul>
        </li> -->
      </ul>
      <div class="queryContainer">
      <form class="d-flex justify-content-end">
        <div class="col-auto">
          <!-- <label for="exampleFormControlSelect1"></label> -->
          <select class="form-select form-select-sm mt-1">
            <option selected>Gene ID</option>
            <option>Uniprot ID</option>
            <option>Keyword</option>
          </select>
        </div>
        <input class="form-control form-control-sm me-1 mt-1 mb-1" type="search" placeholder="Search" aria-label="Search">
        <button class="btn btn-primary btn-sm mt-1 mb-1" type="submit">Search</button>
        <!-- <button class="btn btn-default" type="submit">
          <i class="glyphicon glyphicon-search"></i>
        </button> -->
      </form>
      </div>
    </div>
  </div>
</nav>
  <!-- <div class="col-md-12">
    <div class="input-group">
      <input type="text" class="form-control" placeholder="Search for...">
      <span class="input-group-btn">
        <button class="btn btn-primary" type="button">
          <span class="glyphicon glyphicon-search" aria-hidden="true"></span>
        </button>
      </span>
    </div>
  </div> -->
  <div class="container-fluid">
    <div class="row">
      <div class="col-md-2 border border-dark">2</div>
      <div class="col-md-2 border border-dark">2</div>
      <div class="col-md-2 border border-dark">2</div>
      <div class="col-md-6 border border-dark">6</div>
    </div>
  </div>

  <div class="table-responsive">
    <table class="table table-striped table-hover">
      <thead>
        <tr>
          <th class="count">#</th>
          <th class="gene_id">Gene ID</th>
          <th class="enzyme_name">Enzyme Name</th>
          <th class="uniprot_id">Uniprot ID</th>
          <th class="num_vus"># Missense VUS</th>
          <th class="cadd_score">Highest CADD score</th>
          <th class="EC_number">EC #</th>
          <th class="alphafold">AlphaFold Protein Structure</th>
        </tr>
      </thead>
    </table>
  </div>

  <footer>This is footer</footer>
  <!-- <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
  <script src="js/bootstrap.min.js"></script>
  <link href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" rel="stylesheet" type="text/css" />
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script> -->
  <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"
    integrity="sha384-ka7Sk0Gln4gmtz2MlQnikT1wXgYsOg+OMhuP+IlRH9sENBO0LRn5q+8nbTov4+1p" crossorigin="anonymous"></script>
</body>

</html>
