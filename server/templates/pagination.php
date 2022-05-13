<?php
$first_item = ($current_page-1) * $num_per_page + 1;
$last_item = ($current_page * $num_per_page) < $num_results ? $current_page * $num_per_page: $num_results;
?>

<div class="contents_header">
  <?php
  if (strpos($_SERVER['SCRIPT_NAME'], 'gene.php') !== false){
    include('templates/description_gene.html');
  } elseif (strpos($_SERVER['SCRIPT_NAME'], 'variant.php') !== false) {
   include('templates/description_variant.html');
  } elseif (strpos($_SERVER['SCRIPT_NAME'], 'gene_ec.php') !== false){
    include('templates/description_gene.html');
  }?>
  <p class="num_items">Items: <span class="badge rounded-pill bg-light text-dark"><?php echo escape($first_item) ?></span> to <span class="badge rounded-pill bg-light text-dark"><?php echo escape($last_item); ?></span> of <span class="badge rounded-pill bg-secondary"><?php echo escape($num_results) ?></span></h2>
  <div class="pageforms">
    <form action="" method="Post" class="pageform">
      <button class="page_button btn btn-outline-primary btn-xs" type="submit" name="page" value="1" <?php if($current_page==1){?> disabled="disabled" <?php } ?>>&lt&lt First</button>
      <button class="page_button btn btn-outline-primary btn-xs" type="submit" name="page" value=<?php echo ($current_page-1); ?> <?php if($current_page==1){?>disabled="disabled" <?php }; ?> >&lt Prev</button>
    </form>
    <!-- <form action="" method="Post" class="pageform" onsubmit="validateInput(<?= $total_page ?>)"> -->
    <form action="" method="Post" class="pageform" onsubmit="return validateInput(<?= $total_page ?>)">
      <span>Page </span><input type="text" id="page_number" name="page" value="<?php echo $current_page; ?>" />/ <span id="total_page"><?php echo $total_page ?></span>
      <input class="btn btn-primary btn-xs" type="submit" id="jump_button" value="Go"/>
    </form>
    <form action="" method="Post" class="pageform">
      <button class="page_button btn btn-outline-primary btn-xs" type="submit" name="page" value=<?php echo ($current_page+1); ?> <?php if($current_page==$total_page){?>disabled="disabled" <?php }; ?> >&gt Next</button>
      <button class="page_button btn btn-outline-primary btn-xs" type="submit" name="page" value=<?php echo $total_page ?> <?php if($current_page==$total_page){?>disabled="disabled" <?php } ?>>&gt&gt Last</button>
    </from>
  </div>
</div>

<script>
    function validateInput(stop){
        let x;
        x = document.getElementById("view_page").value;

        if(isNaN(x) || x < 1 || x > stop){
            alert("Enter the number within the page range");
            return false;
        }

        return true;
    }
</script>
