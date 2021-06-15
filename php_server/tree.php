<?php require "templates/header.php"; ?>
<div class="content">

<h2>Browse by EC number</h2>

<?php 
try{
  require "config.php";
  require "common.php";
  $connection = new PDO($dsn, $username, $password, $options);
  // $sql = "SELECT e.ec_number, e.description, e.class, COUNT(m.mutation_id) AS num_vus
  //         FROM Enzyme_class AS e
  //         LEFT JOIN Gene AS g
  //         ON e.ec_number = g.ec_number
  //         LEFT JOIN Mutation AS m
  //         ON g.gene_id = m.gene_id
  //         GROUP BY ec_number, description, class";
  $sql = "SELECT e.ec_number, e.ec_1, e.ec_2, e.ec_3, e.ec_4, description, class, COUNT(m.mutation_id) AS num_vus
          FROM Enzyme_class AS e
          LEFT JOIN Gene AS g
          ON (class = 1 AND e.ec_1 = g.ec_1)
             OR (class = 2 AND e.ec_1 = g.ec_1 AND e.ec_2 = g.ec_2)
             OR (class = 3 AND e.ec_1 = g.ec_1 AND e.ec_2 = g.ec_2 AND e.ec_3 = g.ec_3)
             OR (class = 4 AND e.ec_1 = g.ec_1 AND e.ec_2 = g.ec_2 AND e.ec_3 = g.ec_3 AND e.ec_4 = g.ec_4)
          LEFT JOIN Mutation AS m USING(gene_id)
          GROUP BY ec_number, e.ec_1, e.ec_2, e.ec_3, e.ec_4, description, class
          ORDER BY e.ec_1 ASC, e.ec_2 ASC, e.ec_3 ASC, e.ec_4 ASC";
  $statement = $connection->prepare($sql);
  $statement->execute();
  $result = $statement->fetchAll(); 
}catch(PDOException $error) {
  echo $sql . "<br>" . $error->getMessage();
}
?>

<?php
$ct_class = 1;
if ($result && $statement->rowCount() > 0) { ?>

<ul id="myUL">

<?php foreach ($result as $row) { 
  while($ct_class > $row["class"]){?>
       </ul> 
       </li>
    <?php $ct_class--; 
  }?>
  <?php if($row["class"] == 1){?>
  <li><?php
    if($row["num_vus"] > 0){?>
    <span class="caret"><a href="gene_ec.php?ec=<?php echo($row["ec_1"]) ?>"><?php echo escape($row["ec_number"]." "); echo($row["description"]);?><?php echo escape(" (".$row["num_vus"].")");?></a>
    </span>
  <?php } else {?>
    <span class="caret"><?php echo escape($row["ec_number"]." "); echo($row["description"]);?><?php echo escape(" (".$row["num_vus"].")");?></span>
  <?php } ?>

  <?php } elseif($row["class"] == 2){
    if($ct_class == 1){?> 
    <ul class="nested">
    <?php $ct_class = 2;}?>
    <li><?php
    if($row["num_vus"] > 0){?>
    <span class="caret"><a href="gene_ec.php?ec=<?php echo($row["ec_1"].".".$row["ec_2"]) ?>"><?php echo escape($row["ec_number"]." "); echo($row["description"]);?><?php echo escape(" (".$row["num_vus"].")");?></a>
    </span>
  <?php } else {?>
    <span class="caret"><?php echo escape($row["ec_number"]." "); echo($row["description"]);?><?php echo escape(" (".$row["num_vus"].")");?></span>
  <?php } ?>
  
  <?php } elseif($row["class"] == 3){
    if($ct_class == 2){?> 
    <ul class="nested">
    <?php $ct_class = 3;}?>
    <li><?php
    if($row["num_vus"] > 0){?>
    <span class="caret">
    <a href="gene_ec.php?ec=<?php echo($row["ec_1"].".".$row["ec_2"].".".$row["ec_3"]) ?>"><?php echo escape($row["ec_number"]." "); echo($row["description"]);?><?php echo escape(" (".$row["num_vus"].")");?></a>
    </span>
  <?php } else {?>
    <span class="caret"><?php echo escape($row["ec_number"]." "); echo($row["description"]);?><?php echo escape(" (".$row["num_vus"].")");?></span>
  <?php } ?>

  <?php } elseif($row["class"] == 4){
    if($ct_class == 3){?> 
    <ul class="nested">
    <?php $ct_class = 4;}?>
    <li><?php
      if($row["num_vus"] > 0){?> 
      <a href="gene_ec.php?ec=<?php echo($row["ec_1"].".".$row["ec_2"].".".$row["ec_3"].".".$row["ec_4"]) ?>"><?php echo escape($row["ec_number"]." "); echo($row["description"]);?></a> <?php echo escape(" (".$row["num_vus"].")");?>
      <?php } else{ ?>
      <?php echo escape($row["ec_number"]." "); echo($row["description"]); echo escape(" (0)")?>
      <?php } ?>
    </li>
  <?php } ?>

<?php } ?>
<?php while($ct_class > $row["class"]){?>
       </ul> 
       </li>
  <?php $ct_class--; 
}?>


<?php } else { ?>
    <p>> No results are available.</p>
  <?php } ?>


<script>
var toggler = document.getElementsByClassName("caret");
// toggler = window.getComputedStyle(toggler, '::before')
var i;

for (i = 0; i < toggler.length; i++) {
  toggler[i].addEventListener("click", function() {
    this.parentElement.querySelector(".nested").classList.toggle("active");
    this.classList.toggle("caret-down");
  });
}
</script>

</div>

<?php require "templates/footer.php"; ?>