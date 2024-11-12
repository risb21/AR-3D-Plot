using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

public class ImageReader : MonoBehaviour {
    Texture2D img;
    void Start() {
        img = Resources.Load<Texture2D>("hi");

        if (img == null) {
            Debug.Log("Image couldn't be loaded");
            return;
        }
        
        Debug.Log($"Image loaded\nRandomPixel: {img.GetPixel(3, 3)}");

    }

    void Update() {
        
    }
}
