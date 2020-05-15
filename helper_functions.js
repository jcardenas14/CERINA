function detect_click(el) {
  Shiny.setInputValue('clicked_text', el.innerText, {priority: "event"});
}

function detect_click_mirna(el) {
  Shiny.setInputValue('clicked_mirna_text', el.innerText, {priority: "event"});
}

function detect_click_function(el) {
  Shiny.setInputValue('clicked_text_function', el.innerText, {priority: "event"});
}

shinyjs.init = function(){
  $('#enrichmentWorkflow li a[data-value="func_enrich_val"]').hide();
}