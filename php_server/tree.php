<!DOCTYPE html>
<html>
<head>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
ul, #myUL {
  list-style-type: none;
}

#myUL {
  margin: 0;
  padding: 0;
}

.caret {
  cursor: pointer;
  -webkit-user-select: none; /* Safari 3.1+ */
  -moz-user-select: none; /* Firefox 2+ */
  -ms-user-select: none; /* IE 10+ */
  user-select: none;
}

.caret::before {
  content: "\25B6";
  color: black;
  display: inline-block;
  margin-right: 6px;
}

.caret-down::before {
  -ms-transform: rotate(90deg); /* IE 9 */
  -webkit-transform: rotate(90deg); /* Safari */'
  transform: rotate(90deg);  
}

.nested {
  display: none;
}

.active {
  display: block;
}
</style>
</head>
<body>

<h2>Browse by EC number</h2>
<p>A tree view represents a hierarchical view of information, where each item can have a number of subitems.</p>
<p>Click on the arrow(s) to open or close the tree branches.</p>

<!-- <ul id="myUL">
  <li><span class="caret">1 Oxidoreductases</span>
    <ul class="nested">
      <li><span class="caret">1.1 Acting on the CH-OH group of donors</span>
        <ul class="nested">
          <li><span class="caret">1.1.1 With NAD+ or NADP+ as acceptor</span>
            <ul class="nested">
              <li>1.1.1.1 alcohol dehydrogenase</li>
              <li>1.1.1.2 ...</li>
              <li>1.1.1.3 ...</li>
              <li>1.1.1.4 ...</li>
            </ul>
          </li>
          <li><span class="caret">1.1.2 ...</span>
            <ul class="nested">
              <li>1.1.2.1 ...</li>
              <li>1.1.2.2 ...</li>
              <li>1.1.2.3 ...</li>
              <li>1.1.2.4 ...</li>
            </ul>
          </li>
          <li><span class="caret">1.1.3 ...</span>
            <ul class="nested">
              <li>1.1.3.1 ...</li>
              <li>1.1.3.2 ...</li>
              <li>1.1.3.3 ...</li>
              <li>1.1.3.4 ...</li>
            </ul>
          </li>
        </ul>
      </li>  
      <li><span class="caret">1.2 Acting on ...</span>
        <ul class="nested">
          <li><span class="caret">1.2.1 ...</span>
            <ul class="nested">
              <li>1.2.1.1 ...</li>
              <li>1.2.1.2 ...</li>
              <li>1.2.1.3 ...</li>
              <li>1.2.1.4 ...</li>
            </ul>
          </li>
        </ul>
      </li>  
    </ul>
  </li>
  <li><span class="caret">2 Transferases</span>
    <ul class="nested">
      <li><span class="caret">2.1 Transferring one-carbon groups</span>
        <ul class="nested">
          <li><span class="caret">2.1.1 Methyltransferases</span>
          </li>
        </ul>
      </li>  
      <li><span class="caret">2.7 Transferring phosphorus-containing groups</span>
        <ul class="nested">
          <li><span class="caret">2.7.11 Protein-serine/threonine kinases</span>
            <ul class="nested">
              <li>2.7.11.1 non-specific serine/threonine protein kinase</li>
              <li>2.7.11.17 Ca2+/calmodulin-dependent protein kinase</li>
            </ul>
          </li>
        </ul>
      </li>  
    </ul>
  </li>
</ul> -->

<?php 
try{
  require "config.php";
  require "common.php";
  $connection = new PDO($dsn, $username, $password, $options);
  $sql = "SELECT ec_number, description, class 
          FROM Enzyme_class";
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
  <li><span class="caret"><?php echo escape($row["ec_number"]." "); echo($row["description"]);?></span>

  <?php } elseif($row["class"] == 2){
    if($ct_class == 1){?> 
    <ul class="nested">
    <?php $ct_class = 2;}?>
    <li><span class="caret"><?php echo escape($row["ec_number"]." "); echo($row["description"]);?></span>

  <?php } elseif($row["class"] == 3){
    if($ct_class == 2){?> 
    <ul class="nested">
    <?php $ct_class = 3;}?>
    <li><span class="caret"><?php echo escape($row["ec_number"]." "); echo($row["description"]); ?></span>
  
  <?php } elseif($row["class"] == 4){
    if($ct_class == 3){?> 
    <ul class="nested">
    <?php $ct_class = 4;}?>
    <li><?php echo escape($row["ec_number"]." "); echo($row["description"]); ?>
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
var i;

for (i = 0; i < toggler.length; i++) {
  toggler[i].addEventListener("click", function() {
    this.parentElement.querySelector(".nested").classList.toggle("active");
    this.classList.toggle("caret-down");
  });
}
</script>

</body>
</html>

<a href="mutation.php?gene_id=<?php echo $row["gene_id"] ?>"><?php echo escape($row["gene_name_short"]); ?></a><