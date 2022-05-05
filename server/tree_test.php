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

    .contents{
      font-size: 15px;
      width: 90%;
      margin: 0px auto 50px;
    }

    hr{
      color: grey;
    }

    footer{
      width: 90%;
      margin: 50px auto 50px;
      color: grey;
      font-size: 14px;
    }

    .treeview .list-group-item {
        cursor: pointer
    }

    .treeview span.indent {
        margin-left: 10px;
        margin-right: 10px
    }

    .treeview span.icon {
        width: 12px;
        margin-right: 5px
    }

    .treeview .node-disabled {
        color: silver;
        cursor: not-allowed
    }

    .node-treeview1 {
    }

    .node-treeview1:not(.node-disabled):hover {
            background-color: #F5F5F5;
    }

    h2.page-header {
    margin-top: 0px;
    padding-top: 0px;
    line-height: 15px;
    vertical-align: middle;
}

.table-sortable > thead > tr > th {
    cursor: pointer;
    position: relative;
}

.table-sortable > thead > tr > th:after {
    content: ' ';
    position: absolute;
    height: 0;
    width: 0;
    right: 10px;
    top: 16px;
}

.table-sortable > thead > tr > th:after {
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 5px solid #ccc;
    border-bottom: 0px solid transparent;
}

.table-sortable > thead > tr > th:hover:after {
    border-top: 5px solid #888;
}

.table-sortable > thead > tr > th.asc:after {
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 0px solid transparent;
    border-bottom: 5px solid #333;
}
.table-sortable > thead > tr > th.asc:hover:after {
    border-bottom: 5px solid #888;
}

.table-sortable > thead > tr > th.desc:after {
    /* border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 5px solid #333;
    border-bottom: 5px solid transparent; */
}

  </style>
</head>


<body>
<!-- <?php require "templates/header2.php"; ?> -->
<div class="contents">

<div class="container">
    <h2 class="page-header">Bootstrap Table with Sort Indicators</h2>

    <table class="table table-bordered table-sortable">
        <!-- <caption>
            <code>table.table-sortable</code>
        </caption> -->
        <thead>
            <tr>
                <th class="asc">#</th>
                <th class="desc">Field 1</th>
                <th>Field 2</th>
                <th>Field 3</th>
            </tr>
        </thead>
        <tbody>
            <tr><th scope="row">1</th><td>Data 1</td><td>John</td><td>2020-01-24</td></tr>
            <tr><th scope="row">2</th><td>Data 2</td><td>Andrew</td><td>2020-01-22</td></tr>
            <tr><th scope="row">3</th><td>Data 3</td><td >Laura</td><td>2020-01-20</td></tr>
        </tbody>
    </table>
</div>

</div>
</body>

<?php require "templates/footer2.php"; ?>
