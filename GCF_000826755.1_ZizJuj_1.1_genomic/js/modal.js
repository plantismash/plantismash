$(document).ready(function() {
    $("body").append("<div id='ps-modal' class='modal'><div class='modal-window'><span class='modal-close'>&times;</span><div class='modal-content'></div></div></div>");
	var modal = document.getElementById("ps-modal");
	var closeModal = document.getElementsByClassName("modal-close")[0];
	// When the user clicks on <span> (x), close the modal
	closeModal.onclick = function() {
		modal.style.display = "none";
	}
	// When the user clicks anywhere outside of the modal, close it
	window.onclick = function(event) {
		if (event.target == modal) {
			modal.style.display = "none";
		}
	}
});

function openDialog(contentDom) {
	var modal = $("#ps-modal");
	var modalContent = $("#ps-modal .modal-content");
	var content = $("<div>" + contentDom + "</div>");
	$(modalContent).html(content);
	$(modal).css("display", "block");
}