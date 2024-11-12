using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEditor.Rendering;
using UnityEngine;

[RequireComponent(typeof(Rigidbody))]
public class CameraControlller : MonoBehaviour {
    
    public float sensitivity;
    public float slowSpeed, normalSpeed, sprintSpeed;
    float currentSpeed;

    void Start() {
        string a = Directory.GetCurrentDirectory();
        Debug.Log($"Current Dir: {a}");
    }

    void Update() {
        if (Input.GetMouseButton(1)) {
            Cursor.visible = false;
            Cursor.lockState = CursorLockMode.Locked;
            Movement();
            Rotation();
        } else {
            Cursor.visible = true;
            Cursor.lockState = CursorLockMode.None;
        }
    }

    void Rotation() {
        Vector3 mouseInput = new Vector3(-Input.GetAxis("Mouse Y"), Input.GetAxis("Mouse X"), 0);
        transform.Rotate(mouseInput * sensitivity * Time.deltaTime * 50);
        Vector3 eulerRotation = transform.rotation.eulerAngles;
        transform.rotation = Quaternion.Euler(eulerRotation.x, eulerRotation.y, 0);
    }

    void Movement() {
        Vector3 input = new Vector3(Input.GetAxis("Horizontal"), 0f, Input.GetAxis("Vertical"));
        if (Input.GetKey(KeyCode.LeftShift)) {
            currentSpeed = sprintSpeed;
        } else if (Input.GetKey(KeyCode.LeftControl)) {
            currentSpeed = slowSpeed;
        } else {
            currentSpeed = normalSpeed;
        }

        transform.Translate(input * currentSpeed * Time.deltaTime);
    }


}
